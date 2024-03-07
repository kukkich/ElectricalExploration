using InverseTask.DirectTask;
using InverseTask.DirectTask.EquationSystem;
using Microsoft.Extensions.Logging;
using SharpMath;
using SharpMath.FiniteElement._2D;
using SharpMath.FiniteElement.Materials.HarmonicWithoutChi;
using SharpMath.FiniteElement.Materials.Providers;
using SharpMath.Geometry;
using SharpMath.Geometry._2D;
using SharpMath.Matrices;
using SharpMath.Vectors;

namespace InverseTask;

public class FunctionalOptimizer : Method<FunctionalOptimizerConfig>
{
    private readonly IDirectSolver _directSolver;
    private readonly IParameterDirectionSLAESolver _slaeSolver;

    // [i, k, n] -- R derivative
    // by i-th parameter
    // at k-th frequency
    // at j-th source
    private double[][][] _deviatedNumericalFunc = null!;
    // [k, j] contains numeric function values
    // at k-th frequency
    // in j-th measuring point
    private double[][] _numericalFunc;
    private double[][] _nextNumericalFunc;

    private Vector _currentSigma;
    private Vector _nextSigma;

    private ParameterDirectionEquation _equation = null!;
    private Material[] _materials = null!;
    public FunctionalOptimizer(
        FunctionalOptimizerConfig config,
        ILogger logger,
        IDirectSolver directSolver,
        IParameterDirectionSLAESolver slaeSolver
    ) : base(config, logger)
    {
        _directSolver = directSolver;
        _slaeSolver = slaeSolver;
    }

    public void Solve(Grid<Point, Element> grid, Vector sigmaCurrent)
    {
        Allocate();
        var currentMaterialProvider = new FromArrayMaterialProvider(_materials);

        // direct task
        foreach (var (frequency, frequencyIndex) in Config.Frequencies)
        {
            _numericalFunc[frequencyIndex] = _directSolver.Solve(
                grid,
                frequency,
                currentMaterialProvider,
                Config.MeasuringPoints,
                resultMemory: _numericalFunc[frequencyIndex]
            );
        }
        var prevFunctional = ComputeFunctional(_currentSigma, _numericalFunc);

        // deviated direct tasks
        for (var parameterIndex = 0; parameterIndex < Config.SigmaInitial.Length; parameterIndex++)
        {

            var sigma = _currentSigma[parameterIndex];
            var deviation = sigma / 10;

            foreach (var (frequency, frequencyIndex) in Config.Frequencies)
            {
                var deviatedMaterialProvider = new SigmaDeviatedMaterialProvider(
                    currentMaterialProvider,
                    deviation,
                    parameterIndex + 1 // +1 cause 0-th material is air
                );

                _deviatedNumericalFunc[parameterIndex][frequencyIndex] = _directSolver.Solve(
                    grid,
                    frequency,
                    deviatedMaterialProvider,
                    Config.MeasuringPoints,
                    resultMemory: _deviatedNumericalFunc[parameterIndex][frequencyIndex]
                );
            }
        }

        // Fill equation
        for (var i = 0; i < Config.SigmaInitial.Length; i++)
        {
            _equation.F[i] = 0;
            for (var j = 0; j < Config.SigmaInitial.Length; j++)
            {
                _equation.A[i, j] = 0;
            }

            for (var frequencyIndex = 0; frequencyIndex < Config.Frequencies.Length; frequencyIndex++)
            {
                for (var pointIndex = 0; pointIndex < Config.MeasuringPoints.Length; pointIndex++)
                {
                    var ithParamDiff = -1d * (
                        _deviatedNumericalFunc[i][frequencyIndex][pointIndex] -
                        _numericalFunc[frequencyIndex][pointIndex]
                    ) / (_currentSigma[i] / 10);
                    for (var j = 0; j < Config.SigmaInitial.Length; j++)
                    {

                        var jthParamDiff = -1d * (
                            _deviatedNumericalFunc[j][frequencyIndex][pointIndex] -
                            _numericalFunc[frequencyIndex][pointIndex]
                        ) / (_currentSigma[j] / 10);

                        // TODO implement w coefficient
                        _equation.A[i, j] += ithParamDiff * jthParamDiff;
                    }

                    var measuringDiff = Config.Measurements[frequencyIndex, pointIndex] -
                                        _numericalFunc[frequencyIndex][pointIndex];
                    // TODO implement w coefficient
                    _equation.F[i] = ithParamDiff * measuringDiff;
                }
            }
        }

        var sigmaSeekDirection = _slaeSolver.Solve(_equation);
        _nextSigma = LinAl.LinearCombination(
            _currentSigma, sigmaSeekDirection, 
            1d, Config.Betta, 
            resultMemory: _nextSigma
        );

        // direct task
        foreach (var (frequency, frequencyIndex) in Config.Frequencies)
        {
            _nextNumericalFunc[frequencyIndex] = _directSolver.Solve(
                frequency,
                currentMaterialProvider,
                Config.MeasuringPoints,
                resultMemory: _nextNumericalFunc[frequencyIndex]
            );
        }
        var nextFunctional = ComputeFunctional(_nextSigma, _nextNumericalFunc);

        if (nextFunctional >= prevFunctional)
        {
            Config.Betta /= 2;
        }
    }

    private double ComputeFunctional(Vector sigma, double[][] func)
    {
        var funcPunish = 0d;
        for (var frequencyIndex = 0; frequencyIndex < Config.Frequencies.Length; frequencyIndex++)
        {
            for (var pointIndex = 0; pointIndex < Config.MeasuringPoints.Length; pointIndex++)
            {
                var measuringDiff = Config.Measurements[frequencyIndex, pointIndex] -
                                    func[frequencyIndex][pointIndex];
                // TODO implement w coefficient
                funcPunish += Math.Pow(measuringDiff, 2);
            }
        }

        var sigmaPunish = 0d;
        for (var parameterIndex = 0; parameterIndex < Config.SigmaInitial.Length; parameterIndex++)
        {
            sigmaPunish += Config.Alpha[parameterIndex] * Math.Pow(sigma[parameterIndex] - Config.SigmaInitial[parameterIndex], 2);
        }

        var functional = funcPunish + sigmaPunish;
        Logger.LogInformation(
            "Functional: {functional:E7} funcPunish: {funcPunish:E7} sigmaPunish: {sigmaPunish:E7}",
            functional, funcPunish, sigmaPunish
        );

        return functional;
    }

    private void Allocate(Grid<Point, Element> grid)
    {
        _deviatedNumericalFunc = new double[Config.SigmaInitial.Length][][];
        for (var parameter = 0; parameter < _deviatedNumericalFunc.Length; parameter++)
        {
            _deviatedNumericalFunc[parameter] = new double[Config.Frequencies.Length][];
            for (var frequency = 0; frequency < _deviatedNumericalFunc[0].Length; frequency++)
            {
                _deviatedNumericalFunc[parameter][frequency] = new double[Config.MeasuringPoints.Length];
            }
        }
        _numericalFunc = new double[Config.Frequencies.Length][];
        _nextNumericalFunc = new double[Config.Frequencies.Length][];
        for (var i = 0; i < Config.Frequencies.Length; i++)
        {
            _numericalFunc[i] = new double[Config.MeasuringPoints.Length];
            _nextNumericalFunc[i] = new double[Config.MeasuringPoints.Length];
        }

        _currentSigma = Vector.Create(Config.SigmaInitial.Length, i => Config.SigmaInitial[i]);
        _nextSigma = Vector.Create(Config.SigmaInitial.Length);

        _equation = new ParameterDirectionEquation
        {
            A = new Matrix(new double[Config.SigmaInitial.Length, Config.SigmaInitial.Length]),
            Alpha = Config.Alpha,
            F = Vector.Create(Config.SigmaInitial.Length),
            Solution = Vector.Create(Config.SigmaInitial.Length),
            ParameterCurrent = _currentSigma,
            ParameterInitial = Config.SigmaInitial,
        };

        _materials = new Material[Config.SigmaInitial.Length + 1];
        _materials[0] = new Material(FunctionalOptimizerConfig.Lambda, 0);
        for (var i = 1; i < _materials.Length; i++)
        {
            _materials[i] = new Material(FunctionalOptimizerConfig.Lambda, Config.SigmaInitial[i]);
        }
    }
}

public class FunctionalOptimizerConfig
{
    public Matrix Measurements { get; set; } = null!;
    public Vector MeasuringPoints { get; set; } = null!;
    public Vector Frequencies { get; set; } = null!;
    public Vector SigmaInitial { get; set; } = null!;
    public Vector Alpha { get; set; } = null!;
    public double Betta { get; set; }
    public int MaxIteration { get; set; }
    public double Precision { get; set; }

    public const double Lambda = 1d / (4d * Math.PI * 1e-7);
}
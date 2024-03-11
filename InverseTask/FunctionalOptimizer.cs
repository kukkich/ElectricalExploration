using InverseTask.DirectTask;
using InverseTask.DirectTask.EquationSystem;
using InverseTask.EquationSystem;
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
    private int FixedMaterialsCount => _materials.Length - SigmaInitial.Length;

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

    private Vector MeasuringPoints { get; set; } = null!;
    private Matrix Measurements { get; set; } = null!;
    private Vector Frequencies { get; set; } = null!;
    private Vector SigmaInitial { get; set; } = null!;
    private Vector Alpha { get; set; } = null!;


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

    public void Solve(
        Grid<Point, Element> grid, 
        Vector measuringPoints,
        Matrix measurements,
        Vector frequencies,
        Vector sigmaInitial,
        Vector alpha,
        Material[] fixedMaterials
        )
    {
        MeasuringPoints = measuringPoints;
        Measurements = measurements;
        Frequencies = frequencies;
        SigmaInitial = sigmaInitial;
        Alpha = alpha;

        Allocate(grid, fixedMaterials);
        var currentMaterialProvider = new FromArrayMaterialProvider(_materials);

        // direct task
        foreach (var (frequency, frequencyIndex) in Frequencies)
        {
            _numericalFunc[frequencyIndex] = _directSolver.Solve(
                frequency,
                currentMaterialProvider,
                MeasuringPoints,
                resultMemory: _numericalFunc[frequencyIndex]
            );
        }
        var prevFunctional = ComputeFunctional(_currentSigma, _numericalFunc);

        // deviated direct tasks
        for (var parameterIndex = 0; parameterIndex < SigmaInitial.Length; parameterIndex++)
        {

            var sigma = _currentSigma[parameterIndex];
            var deviation = sigma / 10;

            foreach (var (frequency, frequencyIndex) in Frequencies)
            {
                var deviatedMaterialProvider = new SigmaDeviatedMaterialProvider(
                    currentMaterialProvider,
                    deviation,
                    FixedMaterialsCount + parameterIndex
                );

                _deviatedNumericalFunc[parameterIndex][frequencyIndex] = _directSolver.Solve(
                    frequency,
                    deviatedMaterialProvider,
                    MeasuringPoints,
                    resultMemory: _deviatedNumericalFunc[parameterIndex][frequencyIndex]
                );
            }
        }

        // Fill equation
        for (var i = 0; i < SigmaInitial.Length; i++)
        {
            _equation.F[i] = 0;
            for (var j = 0; j < SigmaInitial.Length; j++)
            {
                _equation.A[i, j] = 0;
            }

            for (var frequencyIndex = 0; frequencyIndex < Frequencies.Length; frequencyIndex++)
            {
                for (var pointIndex = 0; pointIndex < MeasuringPoints.Length; pointIndex++)
                {
                    var ithParamDiff = -1d * (
                        _deviatedNumericalFunc[i][frequencyIndex][pointIndex] -
                        _numericalFunc[frequencyIndex][pointIndex]
                    ) / (_currentSigma[i] / 10);
                    for (var j = 0; j < SigmaInitial.Length; j++)
                    {

                        var jthParamDiff = -1d * (
                            _deviatedNumericalFunc[j][frequencyIndex][pointIndex] -
                            _numericalFunc[frequencyIndex][pointIndex]
                        ) / (_currentSigma[j] / 10);

                        // TODO implement w coefficient
                        _equation.A[i, j] += ithParamDiff * jthParamDiff;
                    }

                    var measuringDiff = Measurements[frequencyIndex, pointIndex] -
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
        foreach (var (frequency, frequencyIndex) in Frequencies)
        {
            _nextNumericalFunc[frequencyIndex] = _directSolver.Solve(
                frequency,
                currentMaterialProvider,
                MeasuringPoints,
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
        for (var frequencyIndex = 0; frequencyIndex < Frequencies.Length; frequencyIndex++)
        {
            for (var pointIndex = 0; pointIndex < MeasuringPoints.Length; pointIndex++)
            {
                var measuringDiff = Measurements[frequencyIndex, pointIndex] -
                                    func[frequencyIndex][pointIndex];
                // TODO implement w coefficient
                funcPunish += Math.Pow(measuringDiff, 2);
            }
        }

        var sigmaPunish = 0d;
        for (var parameterIndex = 0; parameterIndex < SigmaInitial.Length; parameterIndex++)
        {
            sigmaPunish += Alpha[parameterIndex] * Math.Pow(sigma[parameterIndex] - SigmaInitial[parameterIndex], 2);
        }

        var functional = funcPunish + sigmaPunish;
        Logger.LogInformation(
            "Functional: {functional:E7} funcPunish: {funcPunish:E7} sigmaPunish: {sigmaPunish:E7}",
            functional, funcPunish, sigmaPunish
        );

        return functional;
    }

    private void Allocate(Grid<Point, Element> grid, Material[] fixedMaterials)
    {
        _deviatedNumericalFunc = new double[SigmaInitial.Length][][];
        for (var parameter = 0; parameter < _deviatedNumericalFunc.Length; parameter++)
        {
            _deviatedNumericalFunc[parameter] = new double[Frequencies.Length][];
            for (var frequency = 0; frequency < _deviatedNumericalFunc[0].Length; frequency++)
            {
                _deviatedNumericalFunc[parameter][frequency] = new double[MeasuringPoints.Length];
            }
        }
        _numericalFunc = new double[Frequencies.Length][];
        _nextNumericalFunc = new double[Frequencies.Length][];
        for (var i = 0; i < Frequencies.Length; i++)
        {
            _numericalFunc[i] = new double[MeasuringPoints.Length];
            _nextNumericalFunc[i] = new double[MeasuringPoints.Length];
        }

        _currentSigma = Vector.Create(SigmaInitial.Length, i => SigmaInitial[i]);
        _nextSigma = Vector.Create(SigmaInitial.Length);

        _equation = new ParameterDirectionEquation
        {
            A = new Matrix(new double[SigmaInitial.Length, SigmaInitial.Length]),
            Alpha = Alpha,
            F = Vector.Create(SigmaInitial.Length),
            Solution = Vector.Create(SigmaInitial.Length),
            ParameterCurrent = _currentSigma,
            ParameterInitial = SigmaInitial,
        };

        _materials = new Material[fixedMaterials.Length + SigmaInitial.Length];
        for (var i = 0; i < fixedMaterials.Length; i++)
        {
            _materials[i] = fixedMaterials[i];
        }
        for (var i = 0; i < SigmaInitial.Length; i++)
        {
            _materials[fixedMaterials.Length + i] = new Material(FunctionalOptimizerConfig.Lambda, SigmaInitial[i]);
        }

        _directSolver.Allocate(grid);
    }
}

public class FunctionalOptimizerConfig
{
    public double Betta { get; set; }
    public int MaxIteration { get; set; }
    public double Precision { get; set; }

    public const double Lambda = 1d / (4d * Math.PI * 1e-7);
}
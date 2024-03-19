using InverseTask.DirectTask;
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
    public const double DerivativeStepScale = 0.2;
    public const double MaxParameterChangeRatio = 1.5;
    public const double AlphaChangeRation = 1.5;

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
    private FromArrayMaterialProvider _materialProvider;
    private double _currentFunctional;
    private double _functionalChangeRation;

    private Vector MeasuringPoints { get; set; } = null!;
    private Matrix Measurements { get; set; } = null!;
    private Vector Frequencies { get; set; } = null!;
    private Vector SigmaInitial { get; set; } = null!;
    private Vector Alpha { get; set; } = null!;
    private Vector AlphaInitial { get; set; } = null!;

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
        Initialize(grid, measuringPoints, measurements, frequencies, sigmaInitial, alpha, fixedMaterials);

        // direct task
        _currentFunctional = SolveDirectTask(_numericalFunc, SigmaInitial);
        _functionalChangeRation = double.PositiveInfinity;
        Logger.LogDebug("initial sigma {sigma:E15}", _currentSigma);
        for (var i = 1; i <= Config.MaxIteration &&
                        _currentFunctional > Config.Precision &&
                        Math.Abs(_functionalChangeRation - 1) > 1e-4
                        ; i++)
        {
            AlphaInitial.CopyTo(Alpha);
            Iterate();
            Logger.LogDebug("{iteration}-th successful step. Functional ratio {functionalRation:E7}", i, _functionalChangeRation);
        }
    }

    private void Initialize(Grid<Point, Element> grid, Vector measuringPoints, Matrix measurements, Vector frequencies, Vector sigmaInitial,
        Vector alpha, Material[] fixedMaterials)
    {
        MeasuringPoints = measuringPoints;
        Measurements = measurements;
        Frequencies = frequencies;
        SigmaInitial = sigmaInitial;
        AlphaInitial = alpha;
        Alpha = AlphaInitial.Copy();
        Allocate(grid, fixedMaterials);
        _materialProvider = new FromArrayMaterialProvider(_materials);
    }

    private void Iterate()
    {
        // deviated direct tasks
        SolveDeviatedParameterTasks();
        
        var anyParameterConstraintViolated = true;
        while (anyParameterConstraintViolated)
        {
            var equation = FillParameterDirectionSLAE();

            var sigmaSeekDirection = _slaeSolver.Solve(equation);
            _nextSigma = LinAl.LinearCombination(
                _currentSigma, sigmaSeekDirection,
                1d, Config.Betta,
                resultMemory: _nextSigma
            );

            Logger.LogDebug("deltaSigma: {delta:E7}", sigmaSeekDirection);

            

            anyParameterConstraintViolated = CheckParametersConstraintsAndCorrectAlpha(_currentSigma, _nextSigma);

            #region Betta
            // Todo implement betta changing
            //if (nextFunctional >= _currentFunctional)
            //{
            //    Config.Betta /= 2;
            //    Logger.LogInformation("betta decreased: {betta:E3}", Config.Betta);

            //    // rollback material changes
            //    ChangeMaterial(_currentSigma);
            //}
            #endregion
        }
        // direct task with new material
        var nextFunctional = SolveDirectTask(_nextNumericalFunc, _nextSigma);
        _functionalChangeRation = _currentFunctional / nextFunctional;
        _currentFunctional = nextFunctional;
        Logger.LogInformation("next sigma {sigma:E15}", _nextSigma);

        (_numericalFunc, _nextNumericalFunc) = (_nextNumericalFunc, _numericalFunc);
        (_currentSigma, _nextSigma) = (_nextSigma, _currentSigma);
    }

    private bool CheckParametersConstraintsAndCorrectAlpha(Vector previousParameters, Vector nextParameters)
    {
        var anyConstraintViolated = false;

        for (var parameterIndex = 0; parameterIndex < previousParameters.Length; parameterIndex++)
        {
            var previousParameter = previousParameters[parameterIndex];
            var nextParameter = nextParameters[parameterIndex];
            var changeRatio = previousParameter / nextParameter;

            if (!(double.Max(changeRatio, 1d / changeRatio) > MaxParameterChangeRatio))
            {
                continue;
            }

            anyConstraintViolated = true;
            Alpha[parameterIndex] *= AlphaChangeRation;
        }

        if (anyConstraintViolated)
        {
            Logger.LogTrace("Sigma constraints violated");
            Logger.LogDebug("New Alpha: {sigma}", Alpha);
        }
        else
        {
            Logger.LogTrace("Sigma constraints passed");
        }

        return anyConstraintViolated;
    }

    private ParameterDirectionEquation FillParameterDirectionSLAE()
    {
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
                    var w = 1d / Math.Abs(Measurements[frequencyIndex, pointIndex]);
                    var wSquare = w * w;

                    var ithParamDiff = -1d * (
                        _deviatedNumericalFunc[i][frequencyIndex][pointIndex] -
                        _numericalFunc[frequencyIndex][pointIndex]
                    ) / GetDerivativeStep(_currentSigma[i]);
                    for (var j = 0; j < SigmaInitial.Length; j++)
                    {
                        var jthParamDiff = -1d * (
                            _deviatedNumericalFunc[j][frequencyIndex][pointIndex] -
                            _numericalFunc[frequencyIndex][pointIndex]
                        ) / GetDerivativeStep(_currentSigma[j]);

                        _equation.A[i, j] += wSquare * ithParamDiff * jthParamDiff;
                    }

                    var measuringDiff = Measurements[frequencyIndex, pointIndex] -
                                        _numericalFunc[frequencyIndex][pointIndex];
                    _equation.F[i] -= wSquare * ithParamDiff * measuringDiff;
                }
            }
        }

        return _equation;
    }

    private void ChangeMaterial(Vector parametersSource)
    {
        for (var sigmaIndex = 0; sigmaIndex < parametersSource.Length; sigmaIndex++)
        {
            _materials[FixedMaterialsCount + sigmaIndex].Sigma = parametersSource[sigmaIndex];
        }
    }

    private void SolveDeviatedParameterTasks()
    {
        for (var parameterIndex = 0; parameterIndex < SigmaInitial.Length; parameterIndex++)
        {
            // TODO решение при отклонении 2-го (последнего) параметра не меняется
            var sigma = _currentSigma[parameterIndex];
            var deviation = GetDerivativeStep(sigma);
            foreach (var (frequency, frequencyIndex) in Frequencies)
            {
                var deviatedMaterialProvider = new SigmaDeviatedMaterialProvider(
                    _materialProvider,
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
    }

    private double SolveDirectTask(double[][] resultMemory, Vector parameter)
    {
        ChangeMaterial(parameter);

        foreach (var (frequency, frequencyIndex) in Frequencies)
        {
            resultMemory[frequencyIndex] = _directSolver.Solve(
                frequency,
                _materialProvider,
                MeasuringPoints,
                resultMemory: resultMemory[frequencyIndex]
            );
        }

        return ComputeFunctional(parameter, resultMemory);
    }

    public double GetDerivativeStep(double sigma)
    {
        // return 0.2;
        return sigma * 0.1;
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

                var w = 1d / Math.Abs(Measurements[frequencyIndex, pointIndex]);

                funcPunish += Math.Pow(measuringDiff / w, 2);
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
        Logger.LogInformation("{functional:E7}", functional);
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
        _slaeSolver.Allocate(_equation);
    }
}

public class FunctionalOptimizerConfig
{
    public double Betta { get; set; }
    public int MaxIteration { get; set; }
    public double Precision { get; set; }

    public const double Lambda = 1d / (4d * Math.PI * 1e-7);
}
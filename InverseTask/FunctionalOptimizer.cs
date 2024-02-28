using InverseTask.DirectTask.EquationSystem;
using Microsoft.Extensions.Logging;
using SharpMath;
using SharpMath.FiniteElement._2D;
using SharpMath.Geometry;
using SharpMath.Geometry._2D;
using SharpMath.Matrices;
using SharpMath.Vectors;

namespace InverseTask;

public class FunctionalOptimizer : Method<FunctionalOptimizerConfig>
{
    private readonly IDirectSolver _directSolver;
    private readonly IParameterDirectionSLAESolver _slaeSolver;

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

    public void Solve(Grid<Point, Element> grid)
    {

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
}
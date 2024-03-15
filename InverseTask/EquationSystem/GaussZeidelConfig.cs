using Microsoft.Extensions.Logging;
using SharpMath;
using SharpMath.Matrices;
using SharpMath.Vectors;

namespace InverseTask.EquationSystem;

public class GaussZeidelSolver : Method<GaussZeidelConfig>, IAllocationRequired<int>
{
    private Vector _nextSolution;
    private Vector _currentSolution;

    public GaussZeidelSolver(GaussZeidelConfig config, ILogger logger)
        : base(config, logger)
    {
        
    }

    public void Allocate(int dimensionSize)
    {
        _nextSolution = Vector.Create(dimensionSize);
        _currentSolution = Vector.Create(dimensionSize);
    }

    public Vector Solve(MatrixBase matrix, IReadonlyVector<double> rightSide)
    {
        throw new NotImplementedException();
    }
}

public class GaussZeidelConfig
{
    public int MaxIteration { get; set; } = 1000;
    public double Precision { get; set; } = 1e-8;
    public double Relaxation { get; set; } = 1;
}
using SharpMath.Vectors;

namespace InverseTask.EquationSystem;

public class ParameterDirectionSLAESolver : IParameterDirectionSLAESolver
{
    private readonly GaussZeidelSolver _solver;
    private ComputedVector _rightSide;
    private ComputedMatrix _matrix;

    public ParameterDirectionSLAESolver(GaussZeidelSolver solver)
    {
        _solver = solver;
    }

    public Vector Solve(ParameterDirectionEquation equation)
    {
        var solution = _solver.Solve(_matrix, _rightSide);

        return solution;
    }

    public void Allocate(ParameterDirectionEquation param)
    {
        _matrix = new ComputedMatrix(param);
        _rightSide = new ComputedVector(param);
        _solver.Allocate(_rightSide.Length);
    }
}
using SharpMath.Vectors;

namespace InverseTask.DirectTask.EquationSystem;

public interface IParameterDirectionSLAESolver
{
    public Vector Solve(ParameterDirectionEquation equation);
}
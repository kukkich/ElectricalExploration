using InverseTask.DirectTask.EquationSystem;
using SharpMath.Vectors;

namespace InverseTask.EquationSystem;

public interface IParameterDirectionSLAESolver
{
    public Vector Solve(ParameterDirectionEquation equation);
}
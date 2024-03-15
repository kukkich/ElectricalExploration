using SharpMath;
using SharpMath.Vectors;

namespace InverseTask.EquationSystem;

public interface IParameterDirectionSLAESolver : IAllocationRequired<ParameterDirectionEquation>
{
    public Vector Solve(ParameterDirectionEquation equation);
}
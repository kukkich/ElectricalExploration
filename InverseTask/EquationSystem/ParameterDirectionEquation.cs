using SharpMath.Matrices;
using SharpMath.Vectors;

namespace InverseTask.EquationSystem;

public class ParameterDirectionEquation
{
    public Vector Alpha { get; set; } = null!;
    public Matrix A { get; set; } = null!;
    public Vector F { get; set; } = null!;
    public Vector Solution { get; set; } = null!;
    public Vector ParameterInitial { get; set; } = null!;
    public Vector ParameterCurrent { get; set; } = null!;
}
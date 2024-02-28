using SharpMath.Matrices;
using SharpMath.Vectors;

namespace InverseTask.DirectTask.EquationSystem;

public class ParameterDirectionEquation
{
    public Vector Alpha { get; set; }
    public Matrix A { get; set; }
    public Vector F { get; set; }
    public Vector ParameterInitial { get; set; }
    public Vector ParameterCurrent { get; set; }
}
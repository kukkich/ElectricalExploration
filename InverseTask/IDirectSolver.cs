using SharpMath.FiniteElement._2D;
using SharpMath.Geometry;
using SharpMath.Geometry._2D;

namespace InverseTask;

public interface IDirectSolver
{
    public void Solve(Grid<Point, Element> grid);
}
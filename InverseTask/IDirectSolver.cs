using SharpMath.FiniteElement._2D;
using SharpMath.FiniteElement.Core.Assembling.Params;
using SharpMath.FiniteElement.Materials.HarmonicWithoutChi;
using SharpMath.Geometry;
using SharpMath.Geometry._2D;
using SharpMath.Vectors;

namespace InverseTask;

public interface IDirectSolver
{
    public Vector Solve(
        Grid<Point, Element> grid,
        double frequency,
        IMaterialProvider<Material> materialProvider,
        Vector measuringPoints,
        Vector resultMemory
    );
}
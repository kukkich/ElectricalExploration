using SharpMath;
using SharpMath.FiniteElement._2D;
using SharpMath.FiniteElement.Core.Assembling.Params;
using SharpMath.FiniteElement.Materials.HarmonicWithoutChi;
using SharpMath.Geometry;
using SharpMath.Geometry._2D;
using SharpMath.Vectors;

namespace InverseTask.DirectTask;

public interface IDirectSolver : IAllocationRequired<Grid<Point, Element>>
{
    public double[] Solve(
        double frequency,
        IMaterialProvider<Material> materialProvider,
        Vector measuringPoints,
        double[] resultMemory
    );
}
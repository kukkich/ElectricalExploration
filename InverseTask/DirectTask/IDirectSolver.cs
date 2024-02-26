using SharpMath.FiniteElement.Core.Assembling.Params;
using SharpMath.FiniteElement.Materials.HarmonicWithoutChi;
using SharpMath.Vectors;

namespace InverseTask.DirectTask;

public interface IDirectSolver
{
    public Vector Solve(IMaterialProvider<Material> materialProvider, double frequency);
}
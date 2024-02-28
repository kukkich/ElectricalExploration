using SharpMath.FiniteElement.Core.Assembling.Params;
using SharpMath.FiniteElement.Materials.HarmonicWithoutChi;

namespace InverseTask.DirectTask;

public class SigmaDeviatedMaterialProvider : IMaterialProvider<Material>
{
    private readonly IMaterialProvider<Material> _defaultProvider;
    private readonly double _deviation;
    private readonly int _deviatedMaterialId;

    public SigmaDeviatedMaterialProvider(
        IMaterialProvider<Material> defaultProvider,
        double deviation,
        int deviatedMaterialId)
    {
        _defaultProvider = defaultProvider;
        _deviation = deviation;
        _deviatedMaterialId = deviatedMaterialId;
    }

    public Material GetById(int materialId)
    {
        if (materialId != _deviatedMaterialId)
        {
            return _defaultProvider.GetById(materialId);
        }

        var defaultMaterial = _defaultProvider.GetById(materialId);
        var deviatedMaterial = defaultMaterial with {Sigma = defaultMaterial.Sigma + _deviation};

        return deviatedMaterial;
    }
}
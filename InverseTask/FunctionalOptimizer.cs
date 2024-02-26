using Microsoft.Extensions.Logging;
using SharpMath;

namespace InverseTask;

public class FunctionalOptimizer : Method<FunctionalOptimizerConfig>
{
    public FunctionalOptimizer(
        FunctionalOptimizerConfig config, 
        ILogger logger
        ) : base(config, logger)
    {
    }
}

public class FunctionalOptimizerConfig
{
    public double Betta { get; set; }
    public double[] Alpha { get; set; }
    public int MaxIteration { get; set; }
}
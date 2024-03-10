using InverseTask;
using InverseTask.DirectTask;
using Microsoft.Extensions.Configuration;
using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Logging.Abstractions;
using Serilog;
using SharpMath.EquationsSystem.Solver;
using SharpMath.FiniteElement.Materials.HarmonicWithoutChi;
using SharpMath.FiniteElement.Materials.MaterialSetter.Areas;
using SharpMath.FiniteElement.Materials.Providers;
using SharpMath.Geometry._2D;
using SharpMath.Geometry.Splitting;
using SharpMath.Matrices;
using SharpMath.Vectors;

void ConfigureServices(IServiceCollection services)
{
    IConfiguration configuration = new ConfigurationBuilder()
        .SetBasePath(Directory.GetCurrentDirectory())
        .AddJsonFile("appsettings.json", optional: false, reloadOnChange: true)
        .Build();

    services.AddScoped<LocalOptimalSchemeConfig>(provider =>
    {
        provider.GetService<IConfiguration>();
        var losConfig = configuration
            .GetSection("App")
            .GetSection("LOS")
            .Get<LocalOptimalSchemeConfig>();

        return losConfig;
    });

    services.AddSingleton(configuration);

    Log.Logger = new LoggerConfiguration()
        .ReadFrom.Configuration(configuration)
        .CreateLogger();

    services.AddLogging(loggingBuilder =>
        loggingBuilder.AddSerilog(dispose: true));
}


var air = new RectSection(
    new Rectangle(
        0, 0,
        10_000, 5_000
    ),
    materialId: 0
);
var ground = new RectSection(
    new Rectangle(
        0, -200_000,
        10_000, 199_000
    ),
    materialId: 1
);
var conductor = new RectSection(
    new Rectangle(
        0, -1_000,
        10_000, 1_000
    ),
    materialId: 2
);
var areas = new AreasMaterialSetterFactory(
    [conductor, ground, air], 
    defaultMaterialIdId: 1
);
const int nestingDegree = 1;

var grid = new GridBuilder()
    .SetXAxis(new AxisSplitParameter(
        [0, 10_000], 
        new UniformSplitter(50 * nestingDegree)
    ))
    .SetYAxis(new AxisSplitParameter(
        [-200_000, -10_000, -1_000, 0, 5_000],
        [
            new ProportionalSplitter(50 * nestingDegree, 1d/1.05),
            new ProportionalSplitter(50 * nestingDegree, 1d/1.05),
            new ProportionalSplitter(100 * nestingDegree, 1d/1.01),
            new ProportionalSplitter(50 * nestingDegree, 1.01)
        ])
    )
    .SetMaterialSetterFactory(areas)
    .Build();

var solver = new DirectSolver(new LocalOptimalSchemeConfig());
solver.Allocate(grid);

const double frequency = 1e2;
const double lambda = 1d / DirectSolver.MagneticConstant;

var materialProvider = new FromArrayMaterialProvider([
    new Material(lambda, 0),
    new Material(lambda, 1e-3),
    new Material(lambda, 1)
]);
var measuringPoints = new Vector(0d, 200, 400, 1000, 1500, 2000, 3300, 3700, 3950, 4000);
var result = new double[measuringPoints.Length];
solver.Solve(frequency, materialProvider, measuringPoints, result);

var optimizer = new FunctionalOptimizer(
    new FunctionalOptimizerConfig {Betta = 1d, MaxIteration = 10, Precision = 1e-4},
    NullLogger.Instance,
    solver,
    null
    );

var measurements = new Matrix(new double[1,result.Length]);
for (var i = 0; i < result.Length; i++)
{
    measurements[0, i] = result[i];
}

optimizer.Solve(
    grid,
    measuringPoints,
    measurements,
    new Vector(frequency),
    new Vector(1e-3, 1),
    new Vector(1.5, 1.5)
);
Console.WriteLine("Hello, World!");
Console.WriteLine(result[0]);
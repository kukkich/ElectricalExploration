using InverseTask;
using InverseTask.DirectTask;
using InverseTask.EquationSystem;
using Microsoft.Extensions.Configuration;
using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Logging;
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

        return losConfig!;
    });
    services.AddScoped<GaussZeidelConfig>(provider =>
    {
        provider.GetService<IConfiguration>();
        var gaussZeidelConfig = configuration
            .GetSection("App")
            .GetSection("GaussZeidel")
            .Get<GaussZeidelConfig>();

        return gaussZeidelConfig!;
    });
    services.AddScoped<FunctionalOptimizerConfig>(provider =>
    {
        provider.GetService<IConfiguration>();
        var functionalOptimizerConfig = configuration
            .GetSection("App")
            .GetSection("FunctionalOptimizer")
            .Get<FunctionalOptimizerConfig>();

        return functionalOptimizerConfig!;
    });
    services.AddSingleton(configuration);
    
    services.AddScoped<GaussZeidelSolver>();
    services.AddScoped<ParameterDirectionSLAESolver>();
    
    Log.Logger = new LoggerConfiguration()
        .ReadFrom.Configuration(configuration)
        .CreateLogger();
    services.AddLogging(loggingBuilder =>
        loggingBuilder.AddSerilog(dispose: true));
}

void RunSimpleTest()
{
    var services = new ServiceCollection();
    ConfigureServices(services);
    var provider = services.BuildServiceProvider();
    var xAreaSize = 10_000d;
    var xMax = 10_000d;
    var air = new RectSection(
        new Rectangle(
            0, 0,
            xMax, 5_000
        ),
        materialId: 0
    );
    var ground = new RectSection(
        new Rectangle(
            0, -200_000,
            xMax, 199_000
        ),
        materialId: 1
    );
    var conductor = new RectSection(
        new Rectangle(
            0, -1_000,
            xMax, 1_000
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
            [0, xMax],
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
        new Material(lambda, 1e0)
    ]);
    var measuringPoints = new Vector(4000);
    var result = new double[measuringPoints.Length];
    solver.Solve(frequency, materialProvider, measuringPoints, result);

    var optimizer = new FunctionalOptimizer(
        provider.GetRequiredService<FunctionalOptimizerConfig>(),
        provider.GetRequiredService<ILogger<FunctionalOptimizer>>(),
        solver,
        provider.GetRequiredService<ParameterDirectionSLAESolver>()
    );

    var measurements = new Matrix(new double[1, result.Length]);
    for (var i = 0; i < result.Length; i++)
    {
        measurements[0, i] = result[i];
    }

    Console.WriteLine($"R* = {result[0]:E15}");

    optimizer.Solve(
        grid,
        measuringPoints,
        measurements,
        new Vector(frequency),
        sigmaInitial: new Vector(0.1),
        alpha: new Vector(1e-14),
        fixedMaterials: [
            new Material(lambda, 0),
        new Material(lambda, 1e-3)
        ]
    );
    Console.WriteLine(result[0]);
}

void RunFewAreasTest()
{
    var services = new ServiceCollection();
    ConfigureServices(services);
    var provider = services.BuildServiceProvider();
    var xAreaSize = 10_000d;
    var xMax = 30_000d;
    var air = new RectSection(
        new Rectangle(
            0, 0,
            xMax, 5_000
        ),
        materialId: 0
    );
    var ground = new RectSection(
        new Rectangle(
            0, -200_000,
            xMax, 199_000
        ),
        materialId: 1
    );
    List<RectSection> conductors = [
        new RectSection(
            new Rectangle(
                0, -1_000,
                xAreaSize, 1_000
            ),
            materialId: 2
        ),
        new RectSection(
            new Rectangle(
                xAreaSize, -1_000,
                xAreaSize, 1_000
            ),
            materialId: 3
        ),
        new RectSection(
            new Rectangle(
                2*xAreaSize, -1_000,
                xAreaSize, 1_000
            ),
            materialId: 4
        ),
    ];
    
    var areas = new AreasMaterialSetterFactory(
        [conductors[0], conductors[1], conductors[2], ground, air],
        defaultMaterialIdId: 1
    );
    const int nestingDegree = 1;

    var grid = new GridBuilder()
        .SetXAxis(new AxisSplitParameter(
            [0, xAreaSize, 2*xAreaSize, xMax],
            [
                new UniformSplitter(30 * nestingDegree),
                new UniformSplitter(30 * nestingDegree),
                new UniformSplitter(30 * nestingDegree),
            ]
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

    var frequencies = new Vector(new AxisSplitParameter(
            [Math.Pow(10, 0.45d), Math.Pow(10, 2.1d)],
            new UniformSplitter(5)
        )
        .CreateAxis()
        .ToArray());
    const double lambda = 1d / DirectSolver.MagneticConstant;

    var materialProvider = new FromArrayMaterialProvider([
        new Material(lambda, 0),
        new Material(lambda, 1e-3),
        new Material(lambda, 1e-2),
        new Material(lambda, 1e0),
        new Material(lambda, 1e-1)
    ]);
    var measuringPoints = new Vector(new AxisSplitParameter(
            [0, xMax],
            new UniformSplitter(20)
        )
        .CreateAxis()
        .ToArray());
    var result = new double[frequencies.Length][];
    for (var i = 0; i < result.Length; i++)
    {
        result[i] = new double[measuringPoints.Length];
        solver.Solve(frequencies[i], materialProvider, measuringPoints, result[i]);
    }

    var optimizer = new FunctionalOptimizer(
        provider.GetRequiredService<FunctionalOptimizerConfig>(),
        provider.GetRequiredService<ILogger<FunctionalOptimizer>>(),
        solver,
        provider.GetRequiredService<ParameterDirectionSLAESolver>()
    );

    var measurements = new Matrix(new double[frequencies.Length, measuringPoints.Length]);
    for (var i = 0; i < frequencies.Length; i++)
    {
        for (var j = 0; j < measuringPoints.Length; j++)
        {
            measurements[i, j] = result[i][j];
        }
    }
    
    optimizer.Solve(
        grid,
        measuringPoints,
        measurements,
        frequencies,
        sigmaInitial: new Vector(7e-3, 5e-1, 5e-2),
        alpha: Vector.Create(conductors.Count, _ => 1e-14),
        fixedMaterials: [
            new Material(lambda, 0),
            new Material(lambda, 1e-3),
        ]
    );
    Console.WriteLine(result[0]);
}

void RunFewRowsTest()
{
    var services = new ServiceCollection();
    ConfigureServices(services);
    var provider = services.BuildServiceProvider();
    var xAreaSize = 2_000d;
    var xMax = 10_000d;
    var air = new RectSection(
        new Rectangle(
            0, 0,
            xMax, 5_000
        ),
        materialId: 0
    );
    var ground = new RectSection(
        new Rectangle(
            0, -200_000,
            xMax, 199_000
        ),
        materialId: 1
    );
    List<RectSection> firstLayerConductors =
    [
        new RectSection(
            new Rectangle(
                0, -1_000,
                xAreaSize, 1_000
            ),
            materialId: 2
        ),
        new RectSection(
            new Rectangle(
                xAreaSize, -1_000,
                xAreaSize, 1_000
            ),
            materialId: 3
        ),
        new RectSection(
            new Rectangle(
                2 * xAreaSize, -1_000,
                xAreaSize, 1_000
            ),
            materialId: 4
        ),
        new RectSection(
            new Rectangle(
                3 * xAreaSize, -1_000,
                xAreaSize, 1_000
            ),
            materialId: 5
        ),
        new RectSection(
            new Rectangle(
                4 * xAreaSize, -1_000,
                xAreaSize, 1_000
            ),
            materialId: 6
        ),
    ];
    List<RectSection> secondLayerConductors =
    [
        new RectSection(
            new Rectangle(
                0, -2_000,
                xAreaSize, 1_000
            ),
            materialId: 7
        ),
        new RectSection(
            new Rectangle(
                xAreaSize, -2_000,
                xAreaSize, 1_000
            ),
            materialId: 8
        ),
        new RectSection(
            new Rectangle(
                2 * xAreaSize, -2_000,
                xAreaSize, 1_000
            ),
            materialId: 9
        ),
        new RectSection(
            new Rectangle(
                3 * xAreaSize, -2_000,
                xAreaSize, 1_000
            ),
            materialId: 10
        ),
        new RectSection(
            new Rectangle(
                4 * xAreaSize, -2_000,
                xAreaSize, 1_000
            ),
            materialId: 11
        ),
    ];

    var conductors = firstLayerConductors.Union(secondLayerConductors).ToList();
    var areas = conductors.ToList<IMaterialArea<Point>>();
    areas.AddRange([ground, air]);

    var materialSetterFactory = new AreasMaterialSetterFactory(
        areas.ToArray(),
        defaultMaterialIdId: 1
    );
    const int nestingDegree = 1;

    var grid = new GridBuilder()
        .SetXAxis(new AxisSplitParameter(
            [0, xAreaSize, 2 * xAreaSize, 3*xAreaSize, 4*xAreaSize, xMax],
            [
                new UniformSplitter(15 * nestingDegree),
                new UniformSplitter(15 * nestingDegree),
                new UniformSplitter(15 * nestingDegree),
                new UniformSplitter(15 * nestingDegree),
                new UniformSplitter(15 * nestingDegree),
            ]
        ))
        .SetYAxis(new AxisSplitParameter(
            [-200_000, -10_000, -2_000, -1_000, 0, 5_000],
            [
                new ProportionalSplitter(50 * nestingDegree, 1d/1.05),
                new ProportionalSplitter(20 * nestingDegree, 1d/1.05),
                new ProportionalSplitter(30 * nestingDegree, 1d/1.01),
                new ProportionalSplitter(40 * nestingDegree, 1d/1.03),
                new ProportionalSplitter(40 * nestingDegree, 1.03)
            ])
        )
        .SetMaterialSetterFactory(materialSetterFactory)
        .Build();

    var solver = new DirectSolver(new LocalOptimalSchemeConfig());
    solver.Allocate(grid);

    var frequencies = new Vector(new AxisSplitParameter(
            [Math.Pow(10, 0.45d), Math.Pow(10, 2.1d)],
            new ProportionalSplitter(5, 1.1)
        )
        .CreateAxis()
        .ToArray());
    const double lambda = 1d / DirectSolver.MagneticConstant;

    var materialProvider = new FromArrayMaterialProvider([
        new Material(lambda, 0),
        new Material(lambda, 1e-3),

        new Material(lambda, 4e-3),
        new Material(lambda, 3e-2),
        new Material(lambda, 3e-1),
        new Material(lambda, 2e-2),
        new Material(lambda, 1e-3),

        new Material(lambda, 1e-3),
        new Material(lambda, 2e-1),
        new Material(lambda, 1e-0),
        new Material(lambda, 6e-1),
        new Material(lambda, 5e-3),
    ]);
    var measuringPoints = new Vector(new AxisSplitParameter(
            [0, xMax],
            new UniformSplitter(40)
        )
        .CreateAxis()
        .ToArray());
    var result = new double[frequencies.Length][];
    for (var i = 0; i < result.Length; i++)
    {
        result[i] = new double[measuringPoints.Length];
        solver.Solve(frequencies[i], materialProvider, measuringPoints, result[i]);
    }

    var optimizer = new FunctionalOptimizer(
        provider.GetRequiredService<FunctionalOptimizerConfig>(),
        provider.GetRequiredService<ILogger<FunctionalOptimizer>>(),
        solver,
        provider.GetRequiredService<ParameterDirectionSLAESolver>()
    );

    var measurements = new Matrix(new double[frequencies.Length, measuringPoints.Length]);
    for (var i = 0; i < frequencies.Length; i++)
    {
        for (var j = 0; j < measuringPoints.Length; j++)
        {
            measurements[i, j] = result[i][j];
        }
    }
    
    optimizer.Solve(
        grid,
        measuringPoints,
        measurements,
        frequencies,
        sigmaInitial: new Vector(
                1e-3, 7e-2, 8e-1, 4e-1, 5e-3,
                5e-3, 7e-2, 2e-1, 1e-1, 1e-3
            ),
        alpha: Vector.Create(conductors.Count, _ => 1e-14),
        fixedMaterials: [
            new Material(lambda, 0),
            new Material(lambda, 1e-3),
        ]
    );
}

RunFewRowsTest();
using System.Drawing;
using System.Globalization;
using System.Numerics;
using SharpMath;
using SharpMath.EquationsSystem.Solver;
using SharpMath.FiniteElement;
using SharpMath.FiniteElement._2D;
using SharpMath.FiniteElement._2D.Assembling;
using SharpMath.FiniteElement.Core.Assembling;
using SharpMath.FiniteElement.Core.Assembling.Boundary.First;
using SharpMath.FiniteElement.Core.Assembling.Boundary.Second;
using SharpMath.FiniteElement.Core.Harmonic;
using SharpMath.FiniteElement.Materials.MaterialSetter.Areas;
using SharpMath.FiniteElement.Materials.Providers;
using SharpMath.FiniteElement.Providers.Density;
using SharpMath.Geometry;
using SharpMath.Geometry._2D;
using SharpMath.Geometry.Splitting;
using SharpMath.Matrices.Sparse;
using Vector = SharpMath.Vectors.Vector;
using Harmonic2D.Tests.ResultTests;
using Point = SharpMath.Geometry._2D.Point;
using Rectangle = SharpMath.Geometry._2D.Rectangle;
using SharpMath.Matrices.Converters;
using System.Diagnostics;
using SharpMath.EquationsSystem.Preconditions;
using SharpMath.FiniteElement.Core.Harmonic.Solution;

namespace Harmonic2D;

internal class Program
{
    private static object ParallelLock = new();
    static void Main(string[] args)
    {
        Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;
        // RunTest();

        //new ManualTest().Run();

        var testAreas = new Dictionary<string, Rectangle>
        {
            // ["Top\\S"] = new (
            //     1_300d, -700d,
            //     500d, 200d
            // ),
            // ["Top\\M"] = new(
            //     1_300d, -1000d,
            //     500d, 500d
            // ),
            // ["Top\\L"] = new(
            //     1_300d, -1300d,
            //     500d, 800d
            // ),
            // ["Middle\\S"] = new(
            //     1_300d, -2200d,
            //     500d, 200d
            // ),
            // ["Middle\\M"] = new(
            //     1_300d, -2500d,
            //     500d, 500d
            // ),
            // ["Middle\\L"] = new(
            //     1_300d, -2800d,
            //     500d, 800d
            // ),
            // ["Bottom\\S"] = new(
            //     1_300d, -4000d,
            //     500d, 200d
            // ),
            // ["Bottom\\M"] = new(
            //     1_300d, -4300d,
            //     500d, 500d
            // ),
            // ["Bottom\\L"] = new(
            //     1_300d, -4600d,
            //     500d, 800d
            // ),
            //["uniform"] = new (
            //     0d, -5000d,
            //     4000d, 5000d
            // ),
            ["Bottom"] = new(
                9_750d, -2000d,
                500d, 1000d
            ),
            ["Middle"] = new(
                9_750d, -350d,
                500d, 250d
            ),
            ["Top"] = new(
                9_750d, -280d,
                500d, 250d
            ),
        };
        var targetMaterial = 14;
        foreach (var keyValue in testAreas)
        {
            Console.ForegroundColor = ConsoleColor.Green;
            Console.WriteLine($"========================={keyValue.Key}=========================");
            Console.ResetColor();

            var timer = new Stopwatch();
            timer.Start();

            var materialIds = new int[]
            {
                5, 6, 7, 8,
                9, 10, 11, 12,
                13, 14
            };
            foreach (var materialId in materialIds)
            {
                //var testKey = keyValue.Key;
                var areas = new AreasMaterialSetterFactory(new IMaterialArea<Point>[]
                    {
                        
                        // new RectSection(
                        //     new Rectangle(
                        //         500d, -2500d,
                        //         500d, 400d
                        //     ),
                        //     13
                        // ),
                        // new RectSection(
                        //     new Rectangle(
                        //         700d, -4000d,
                        //         600d, 800d
                        //     ),
                        //     13
                        // ),
                        new RectSection(keyValue.Value, materialId),
                        new RectSection(
                            new Rectangle(
                                0d, 0d,
                                40_000d, 5000d
                            ),
                            0
                        ),
                    },
                    defaultMaterialIdId: 1
                );
                RunTestAndWriteResults(materialId, areas, keyValue.Key);
            }
            timer.Stop();
            var ts = timer.Elapsed;
            string elapsedTime = String.Format("{0:00}:{1:00}.{2:00}",
                ts.Minutes, ts.Seconds,
                ts.Milliseconds / 10);

            Console.ForegroundColor = ConsoleColor.Green;
            Console.WriteLine($"{keyValue.Key} Finished {elapsedTime}");
            Console.ResetColor();
        }

        Console.ReadKey();
    }

    private static void RunTestAndWriteResults(int targetMaterial, AreasMaterialSetterFactory areas, string folder)
    {
        int maxDegreeOfParallelism = 1;
        var ws = new[]
        {
             //1e-4, 2e-4,
             //5e-4, 7e-4,
             1e-3, 2e-3,
             5e-3, 7e-3,
             1e-2,
             5e-2,
             1e-1, 5e-1,
             1e-0, 2e-0, 3e-0, 5e-0, 7e-0, 9e-0,
             1e+1, 2e+1, 3e+1, 5e+1, 7e+1, 9e+1,
             1e+2, 2e+2, 3e+2, 5e+2, 7e+2, 9e+2,
             1e+3,
             //2e+3, 5e+3,
             //7e+3,
        };
        for (int i = 0; i < ws.Length; i++)
        {
            ws[i] /= 2 * Math.PI;
        }
        var hs = new double[]
        {
            //9000, 13000, 18000, 
            //30000, 
            //50000,
            200_000
        };
        var xs = new double[]
        {
            10_000,
            //10_250,
            //10_500,
            //10_750,
            //11_000,
            //11_500,
            //12_000,
            //13_000,
            //14_000,
            //14_000,
            //17_000,
            //20_000,
            //30_000,
            //40_000
        };

        foreach (var h in hs)
        {
            var pathBase = "C:\\Users\\vitia\\PycharmProjects\\InverseTasks\\finnalyTetst\\2D Graphs\\Sigma\\"
                           + folder
                           + "\\";
            Console.WriteLine(pathBase);

            var material = LayersTest.Materials;
            var sigma = material[targetMaterial].Sigma;
            var sigmaPathPostFix = sigma.ToString("E2")
                .Replace(".", String.Empty)
                .Replace("00", String.Empty);
            Console.WriteLine(sigmaPathPostFix);

            var options = new ParallelOptions { MaxDegreeOfParallelism = maxDegreeOfParallelism };
            var results = new List<double[]>(xs.Length);
            for (var i = 0; i < xs.Length; i++)
            {
                results.Add(new double[ws.Length]);
            }

            var progress = 0;
            Parallel.For(0, ws.Length, options, i =>
            {
                var testResult = RunResultTest(ws[i], areas, h);
                for (var j = 0; j < xs.Length; j++)
                {
                    var z = testResult.Calculate(new Point(xs[j], 0));
                    var R = z * z * LayersTest.Lambda / ws[i];
                    if (double.IsNaN(R))
                    {
                        throw new Exception("Вырожденное решение");
                    }

                    results[j][i] = R;
                }

                lock (ParallelLock)
                {
                    Console.WriteLine($"[{++progress}/{ws.Length}] {ws[i]:E1} Complete");
                }
            });
            for (int i = 0; i < xs.Length; i++)
            {
                var fullPath = pathBase + $"sigma={sigmaPathPostFix}_" + $"x={xs[i]}.txt";
                using var stream = new StreamWriter(fullPath);
                Console.WriteLine($"Writing...");
                stream.WriteLine(xs[i]);
                WriteValuesString(stream, ws);
                WriteValuesString(stream, results[i]);

                Console.WriteLine($"Written.");
            }
        }
    }

    private static void WriteValuesString(StreamWriter stream, double[] a)
    {
        for (int i = 0; i < a.Length - 1; i++)
        {
            stream.Write($"{a[i]:E15}");
            stream.Write(", ");
        }
        stream.WriteLine($"{a[^1]:E15}");
    }

    private static ImpedanceSolution RunResultTest(double w, AreasMaterialSetterFactory areas, double h)
    {
        const int groundStepsY = 400;
        const int airStepsY = 100;
        const int xSteps = 50;

        var grid = new GridBuilder()
            .SetYAxis(new AxisSplitParameter(new double[]
            {
                -h, 0, 5000
            },
                new ProportionalSplitter(groundStepsY, 1d / 1.05),
                new ProportionalSplitter(airStepsY, 1.1)
            ))
            .SetXAxis(new AxisSplitParameter(new double[]
            {
                0, 9750, 40_000
            },
                new ProportionalSplitter(xSteps, 1d / 1.5),
                new ProportionalSplitter(xSteps, 1.5)
            ))
            .SetMaterialSetterFactory(areas)
            .Build();

        var test = new LayersTest()
            .SetFrequency(w)
            .SetSizes(new Size(2 * xSteps, groundStepsY + airStepsY))
            .SetGrid(grid);

        var context = test.Run();
        var config = new LocalOptimalSchemeConfig
        {
            Eps = 1e-14,
            MaxIterations = 2000
        };

        var profile = MatrixConverter.Convert(context.Equation.Matrix);
        var equationProfile = new Equation<ProfileMatrix>(profile, context.Equation.Solution, context.Equation.RightSide);
        var profileSolver = new LUProfile();
        profileSolver.Solve(equationProfile);
        var solution = new ImpedanceSolution(context.Grid, equationProfile.Solution, w);
        // var preconditioner = new LUPreconditioner();
        // var luSparse = new LUSparse(preconditioner);
        // var solver = new LocalOptimalScheme(preconditioner, luSparse, config);
        // solver.Solve(context.Equation);
        // var solution = new ImpedanceSolution(context.Grid, context.Materials, context.Equation.Solution);
        //
        //  var solution = new FiniteElementSolution2DHarmonic(context.Grid, context.Materials, context.Equation.Solution);
        //
        //  var groundNodes = grid.Nodes.Skip(groundIndexes[0]).Take(groundIndexes.Count).ToArray();


        return solution;
    }

    static EquationAssembler CreateAssembler(HarmonicContext<Point, Element, SparseMatrix> context)
    {
        return new EquationAssembler(
            context,
            new HarmonicLocalAssembler(context),
            new Inserter(),
            new GaussExcluderSparse(),
            new SecondBoundaryApplier(context, new Inserter())
            );
    }

    static Grid<Point, Element> CreateGrid()
    {
        GridBuilder gb = new GridBuilder()
            .SetXAxis(new AxisSplitParameter(new double[]
            {
                0, 2
            }, new IIntervalSplitter[]
            {
                new UniformSplitter(2),
            }))
            .SetYAxis(new AxisSplitParameter(new double[]
            {
                1, 3
            }, new IIntervalSplitter[]
            {
                new UniformSplitter(2),
            }));

        return gb.Build();
    }

    static Context<Point, Element, SparseMatrix> CreateContext(Grid<Point, Element> grid)
    {
        IMatrixPortraitBuilder<SparseMatrix, Element> portraitBuilder = new MatrixPortraitBuilder();
        var matrix = portraitBuilder.Build(grid.Elements, grid.Nodes.Length);

        var context = new Context<Point, Element, SparseMatrix>
        {
            Grid = grid,
            Equation = new Equation<SparseMatrix>(
                Matrix: matrix,
                RightSide: Vector.Create(grid.Nodes.Length * 2),
                Solution: Vector.Create(grid.Nodes.Length * 2)
            ),
            DensityFunctionProvider = null,
            Materials = new DefaultMaterialProvider(),
            FirstConditions = null,
            SecondConditions = null,
        };

        context.DensityFunctionProvider = new AnalyticComplexDensity(context, p => new Complex(-1d * p.X + p.Y, p.X + p.Y));
        context.FirstConditions = CrateFirstConditions(context, p => new Complex(p.X + p.Y, p.X - p.Y));
        return context;
    }

    static FirstCondition[] CrateFirstConditions(Context<Point, Element, SparseMatrix> context, Func<Point, Complex> u)
    {
        var conditions = new FirstBoundaryGeneration[]
        {
            new() {ElementIndex = 0, Bounds = new[] {Bound.Bottom, Bound.Left, }},
            new() {ElementIndex = 1, Bounds = new[] {Bound.Bottom, Bound.Right}},
            new() {ElementIndex = 2, Bounds = new[] {Bound.Top, Bound.Left}},
            new() {ElementIndex = 3, Bounds = new[] {Bound.Top, Bound.Right}},
        };

        var result = new List<FirstCondition>();
        foreach (var condition in conditions)
        {
            var elem = context.Grid.Elements[condition.ElementIndex];
            foreach (var bound in condition.Bounds)
            {
                var nodeIndexes = elem.GetBoundNodeIndexes(bound);
                foreach (var nodeIndex in nodeIndexes)
                {
                    var node = context.Grid.Nodes[nodeIndex];
                    result.Add(new FirstCondition(nodeIndex * 2, u(node).Real));
                    result.Add(new FirstCondition(nodeIndex * 2 + 1, u(node).Imaginary));
                }
            }
        }

        return result.ToArray();
    }

    static SecondCondition[] CrateSecondConditions(Context<Point, Element, SparseMatrix> context, Func<Point, Complex> u)
    {

        return new SecondCondition[]
        {
            new(0, Bound.Bottom, new []{-1d, -1d}, ComponentType.Real)
        };
    }
}

internal class FirstBoundaryGeneration
{
    public int ElementIndex { get; set; }
    public Bound[] Bounds { get; set; }
}
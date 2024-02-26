using System.Drawing;
using System.Numerics;
using SharpMath.FiniteElement;
using SharpMath;
using SharpMath.FiniteElement._2D;
using SharpMath.FiniteElement._2D.Assembling;
using SharpMath.FiniteElement.Core.Assembling;
using SharpMath.FiniteElement.Core.Assembling.Boundary.First;
using SharpMath.FiniteElement.Core.Assembling.Boundary.Second;
using SharpMath.FiniteElement.Materials.HarmonicWithoutChi;
using SharpMath.Geometry;
using SharpMath.Matrices.Sparse;
using SharpMath.FiniteElement.Core.Assembling.Params;
using SharpMath.FiniteElement.Materials.Providers;
using SharpMath.FiniteElement.Providers.Density;
using Vector = SharpMath.Vectors.Vector;
using SharpMath.FiniteElement.Core.Harmonic;
using SharpMath.Geometry._2D;
using Point = SharpMath.Geometry._2D.Point;

namespace Harmonic2D.Tests.ResultTests;

public class LayersTest
{
    public const double Lambda = 1d / (4d * Math.PI * 1e-7);

    private Grid<Point, Element> _grid = null!;
    private Size _sizeInElements;
    private double _w;
    private readonly Material[] _material = Materials.ToArray();
    public static Material[] Materials => new Dictionary<int, Material>
    {
        [0] = new() { Lambda = Lambda, Sigma = 0 },
        [1] = new() { Lambda = Lambda, Sigma = 1e-3 },
        [2] = new() { Lambda = Lambda, Sigma = 2e-3 },
        [3] = new() { Lambda = Lambda, Sigma = 5e-3 },
        [4] = new() { Lambda = Lambda, Sigma = 7e-3 },

        [5] = new() { Lambda = Lambda, Sigma = 1e-2 },
        [6] = new() { Lambda = Lambda, Sigma = 2e-2 },
        [7] = new() { Lambda = Lambda, Sigma = 5e-2 },
        [8] = new() { Lambda = Lambda, Sigma = 7e-2 },

        [9] = new() { Lambda = Lambda, Sigma = 1e-1 },
        [10] = new() { Lambda = Lambda, Sigma = 2e-1 },
        [11] = new() { Lambda = Lambda, Sigma = 5e-1 },
        [12] = new() { Lambda = Lambda, Sigma = 7e-1 },

        [13] = new() { Lambda = Lambda, Sigma = 1 },
        [14] = new() { Lambda = Lambda, Sigma = 2 },
        [15] = new() { Lambda = Lambda, Sigma = 5 },
        [16] = new() { Lambda = Lambda, Sigma = 7 },

        [17] = new() { Lambda = Lambda, Sigma = 10 },
    }.Select(kv => kv.Value).ToArray();
    public LayersTest SetGrid(Grid<Point, Element> grid)
    {
        _grid = grid;
        return this;
    }

    public LayersTest SetFrequency(double w)
    {
        _w = w;
        return this;
    }
    public LayersTest SetSizes(Size size)
    {
        _sizeInElements = size;
        return this;
    }

    public Context<Point, Element, SparseMatrix> Run()
    {
        var context = CreateContext();
        var assembler = CreateAssembler(context);
        assembler.BuildEquation(context)
            .ApplySecondConditions(context)
            .ApplyFirstBoundary(context);

        return context;
    }

    private HarmonicContext<Point, Element, SparseMatrix> CreateContext()
    {
        IMatrixPortraitBuilder<SparseMatrix, Element> portraitBuilder = new MatrixPortraitBuilder();
        var matrix = portraitBuilder.Build(_grid.Elements, _grid.Nodes.Length);

        var context = new HarmonicContext<Point, Element, SparseMatrix>
        {
            Grid = _grid,
            Equation = new Equation<SparseMatrix>(
                Matrix: matrix,
                RightSide: Vector.Create(_grid.Nodes.Length * 2),
                Solution: Vector.Create(_grid.Nodes.Length * 2)
            ),
            Materials = new FromArrayMaterialProvider(_material),
            DensityFunctionProvider = null,
            FirstConditions = null,
            SecondConditions = null,
            Frequency = _w
        };

        context.DensityFunctionProvider = new AnalyticComplexDensity(context, p => new (0,0));
        context.FirstConditions = CreateFirst();
        context.SecondConditions = CreateSecond();

        return context;
    }
    private EquationAssembler CreateAssembler(HarmonicContext<Point, Element, SparseMatrix> context)
    {
        var inserter = new Inserter();

        return new EquationAssembler(
            context,
            new HarmonicLocalAssembler(context),
            inserter,
            new GaussExcluderSparse(),
            new SecondBoundaryApplier(context, inserter)
        );
    }
    private static IEnumerable<int> EnumerateBottomElementIndexes(int xSize, int ySize)
    {
        for (int i = 0; i <= xSize - 1; i++)
        {
            yield return i;
        }
    }
    private static IEnumerable<int> EnumerateTopElementIndexes(int xSize, int ySize)
    {
        for (int i = xSize * (ySize - 1); i <= xSize * ySize - 1; i++)
        {
            yield return i;
        }
    }
    private static IEnumerable<int> EnumerateLeftElementIndexes(int xSize, int ySize)
    {
        for (int i = 0; i <= xSize * (ySize - 1); i += xSize)
        {
            yield return i;
        }
    }
    private static IEnumerable<int> EnumerateRightElementIndexes(int xSize, int ySize)
    {
        for (int i = xSize - 1; i <= xSize * ySize - 1; i += xSize)
        {
            yield return i;
        }
    }

    private SecondCondition[] CreateSecond()
    {
        var result = new List<SecondCondition>();

        foreach (var elementIndex in EnumerateLeftElementIndexes(_sizeInElements.Width, _sizeInElements.Height))
        {
            result.Add(new SecondCondition(elementIndex, Bound.Left, new[] { 0d, 0d }, ComponentType.Real));
            result.Add(new SecondCondition(elementIndex, Bound.Left, new[] { 0d, 0d }, ComponentType.Imaginary));
        }
        foreach (var elementIndex in EnumerateRightElementIndexes(_sizeInElements.Width, _sizeInElements.Height))
        {
            result.Add(new SecondCondition(elementIndex, Bound.Right, new[] { 0d, 0d }, ComponentType.Real));
            result.Add(new SecondCondition(elementIndex, Bound.Right, new[] { 0d, 0d }, ComponentType.Imaginary));
        }

        foreach (var elementIndex in EnumerateTopElementIndexes(_sizeInElements.Width, _sizeInElements.Height))
        {
            result.Add(new SecondCondition(elementIndex, Bound.Bottom, new[] { 1d, 1d }, ComponentType.Real));
            result.Add(new SecondCondition(elementIndex, Bound.Bottom, new[] { 0d, 0d }, ComponentType.Imaginary));
        }

        return result.ToArray();
    }

    private FirstCondition[] CreateFirst()
    {
        var result = new FirstCondition[(_sizeInElements.Width + 1) * 2];

        var i = 0;
        foreach (var elementIndex in EnumerateBottomElementIndexes(_sizeInElements.Width, _sizeInElements.Height))
        {
            var element = _grid.Elements[elementIndex];
            var nodeIndexes = element.GetBoundNodeIndexes(Bound.Bottom);

            foreach (var nodeIndex in nodeIndexes)
            {
                result[nodeIndex * 2] = new FirstCondition(nodeIndex * 2, 0);
                result[nodeIndex * 2 + 1] = new FirstCondition(nodeIndex * 2 + 1, 0);
            }
        }

        return result;
    }
}
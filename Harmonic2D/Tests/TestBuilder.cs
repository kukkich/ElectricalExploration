using System.Diagnostics;
using System.Numerics;
using System.Security.Cryptography.X509Certificates;
using SharpMath;
using SharpMath.FiniteElement;
using SharpMath.FiniteElement._2D;
using SharpMath.FiniteElement._2D.Assembling;
using SharpMath.FiniteElement.Core.Assembling;
using SharpMath.FiniteElement.Core.Assembling.Boundary.First;
using SharpMath.FiniteElement.Core.Assembling.Boundary.Second;
using SharpMath.FiniteElement.Core.Assembling.Params;
using SharpMath.FiniteElement.Core.Harmonic;
using SharpMath.FiniteElement.Materials.HarmonicWithoutChi;
using SharpMath.FiniteElement.Providers.Density;
using SharpMath.Geometry;
using SharpMath.Geometry._2D;
using SharpMath.Matrices.Sparse;
using Vector = SharpMath.Vectors.Vector;

namespace Harmonic2D.Tests;

public class TestBuilder
{
    private Func<Point, double> _fSin;
    private Func<Point, double> _fCos;
    private Func<Point, double> _uSin;
    private Func<Point, double> _uCos;

    private IMaterialProvider<Material> _material;
    private Grid<Point, Element> _grid;
    private BoundaryCreation[] _first;
    private BoundaryCreation[] _second;

    private Func<Point, double> _uSinXDerivative;
    private Func<Point, double> _uSinYDerivative;
    private Func<Point, double> _uCosXDerivative;
    private Func<Point, double> _uCosYDerivative;

    public TestBuilder SetFSin(Func<Point, double> fSin)
    {
        _fSin = fSin;
        return this;
    }
    public TestBuilder SetFCos(Func<Point, double> fCos)
    {
        _fCos = fCos;
        return this;
    }
    public TestBuilder SetUSin(Func<Point, double> uSin)
    {
        _uSin = uSin;
        return this;
    }
    public TestBuilder SetUCos(Func<Point, double> uCos)
    {
        _uCos = uCos;
        return this;
    }

    public TestBuilder SetMaterial(IMaterialProvider<Material> materialProvider)
    {
        _material = materialProvider;
        return this;
    }
    public TestBuilder SetGrid(Grid<Point, Element> grid)
    {
        _grid = grid;
        return this;
    }
    public TestBuilder SetFirstBoundary(BoundaryCreation[] boundary)
    {
        _first = boundary;
        return this;
    }
    public TestBuilder SetSecondBoundary(BoundaryCreation[] boundary)
    {
        _second = boundary;
        return this;
    }

    public TestBuilder SetFirstAround(int xSize, int ySize, ComponentType type)
    {
        _first ??= Array.Empty<BoundaryCreation>();
        _first = _first
            .Union(EnumerateBottomBoundary(xSize, ySize, type))
            .Union(EnumerateRightBoundary(xSize, ySize, type))
            .Union(EnumerateTopBoundary(xSize, ySize, type))
            .Union(EnumerateLeftBoundary(xSize, ySize, type))
            .OrderBy(x => x.ElemIndex)
            .ThenBy(x => x.ComponentType)
            .ToArray();

        return this;
    }

    public TestBuilder SetUSinXDerivative(Func<Point, double> uSin)
    {
        _uSinXDerivative = uSin;
        return this;
    }
    public TestBuilder SetUSinYDerivative(Func<Point, double> uSin)
    {
        _uSinYDerivative = uSin;
        return this;
    }
    public TestBuilder SetUCosXDerivative(Func<Point, double> uCos)
    {
        _uCosXDerivative = uCos;
        return this;
    }
    public TestBuilder SetUCosYDerivative(Func<Point, double> uCos)
    {
        _uCosYDerivative = uCos;
        return this;
    }

    public Context<Point, Element, SparseMatrix> Build()
    {
        var context = CreateContext();
        var assembler = CreateAssembler(context);
        assembler.BuildEquation(context)
            .ApplySecondConditions(context)
            .ApplyFirstBoundary(context);
        
        return context;
    }

    private Context<Point, Element, SparseMatrix> CreateContext()
    {
        IMatrixPortraitBuilder<SparseMatrix, Element> portraitBuilder = new MatrixPortraitBuilder();
        var matrix = portraitBuilder.Build(_grid.Elements, _grid.Nodes.Length);

        var context = new Context<Point, Element, SparseMatrix>
        {
            Grid = _grid,
            Equation = new Equation<SparseMatrix>(
                Matrix: matrix,
                RightSide: Vector.Create(_grid.Nodes.Length * 2),
                Solution: Vector.Create(_grid.Nodes.Length * 2)
            ),
            DensityFunctionProvider = null,
            Materials = _material,
            FirstConditions = null,
            SecondConditions = null,
        };

        context.DensityFunctionProvider = new AnalyticComplexDensity(context, p => new Complex(_fSin(p), _fCos(p)));
        context.FirstConditions = CreateFirst();
        context.SecondConditions = CreateSecond();
        return context;
    }
    private EquationAssembler CreateAssembler(Context<Point, Element, SparseMatrix> context)
    {
        var inserter = new Inserter();

        return new EquationAssembler(
            context,
            new LocalAssembler(context),
            inserter,
            new GaussExcluderSparse(),
            new SecondBoundaryApplier(context, inserter)
        );
    }
    private FirstCondition[] CreateFirst()
    {
        var result = new List<FirstCondition>();

        foreach (var condition in _first)
        {
            var elem = _grid.Elements[condition.ElemIndex];
            foreach (var bound in condition.Bounds)
            {
                var indexes = elem.GetBoundNodeIndexes(bound);
                var points = indexes.Select(x => _grid.Nodes[x]);
                var values = points.Select(p =>
                    condition.ComponentType is ComponentType.Real ? _uSin(p) : _uCos(p)
                ).ToArray();

                for (int i = 0; i < indexes.Length; i++)
                {
                    var index = indexes[i] * 2;
                    if (condition.ComponentType is ComponentType.Imaginary)
                        index += 1;

                    result.Add(new FirstCondition(index, values[i]));
                }
            }
        }

        result = result.DistinctBy(x => x.NodeIndex)
            .OrderBy(x => x.NodeIndex)
            .ToList();
        return result.ToArray();
    }
    private SecondCondition[] CreateSecond()
    {
        return (from condition in _second
                let elem = _grid.Elements[condition.ElemIndex]
                from bound in condition.Bounds
                let u = (condition.ComponentType, bound) switch
                {
                    (ComponentType.Real, Bound.Right) => _uSinXDerivative,
                    (ComponentType.Real, Bound.Left) => p => -1d * _uSinXDerivative(p),
                    (ComponentType.Real, Bound.Top) => _uSinYDerivative,
                    (ComponentType.Real, Bound.Bottom) => p => -1d * _uSinYDerivative(p),
                    (ComponentType.Imaginary, Bound.Right) => _uCosXDerivative,
                    (ComponentType.Imaginary, Bound.Left) => p => -1d * _uCosXDerivative(p),
                    (ComponentType.Imaginary, Bound.Top) => _uCosYDerivative,
                    (ComponentType.Imaginary, Bound.Bottom) => p => -1d * _uCosYDerivative(p),
                }
                let indexes = elem.GetBoundNodeIndexes(bound)
                let points = indexes.Select(x => _grid.Nodes[x])
                let values = points.Select(u).ToArray()
                select new SecondCondition(condition.ElemIndex, bound, values, condition.ComponentType)).ToArray();
    }

    private static IEnumerable<BoundaryCreation> EnumerateBottomBoundary(int xSize, int ySize, ComponentType type)
    {
        for (int i = 0; i <= xSize - 1; i++)
        {
            yield return new BoundaryCreation(i, new[] { Bound.Bottom }, type);
        }
    }
    private static IEnumerable<BoundaryCreation> EnumerateTopBoundary(int xSize, int ySize, ComponentType type)
    {
        for (int i = xSize * (ySize - 1); i <= xSize * ySize - 1; i++)
        {
            yield return new BoundaryCreation(i, new[] { Bound.Top }, type);
        }
    }
    private static IEnumerable<BoundaryCreation> EnumerateLeftBoundary(int xSize, int ySize, ComponentType type)
    {
        for (int i = 0; i <= xSize * (ySize - 1); i += xSize)
        {
            yield return new BoundaryCreation(i, new[] { Bound.Left }, type);
        }
    }
    private static IEnumerable<BoundaryCreation> EnumerateRightBoundary(int xSize, int ySize, ComponentType type)
    {
        for (int i = xSize - 1; i <= xSize * ySize - 1; i += xSize)
        {
            yield return new BoundaryCreation(i, new[] { Bound.Right }, type);
        }
    }
}

[DebuggerDisplay("{ElemIndex} {Bounds[0]} {ComponentType}")]
public record struct BoundaryCreation(int ElemIndex, Bound[] Bounds, ComponentType ComponentType);
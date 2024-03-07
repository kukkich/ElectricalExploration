using System.Numerics;
using InverseTask.DirectTask.Result;
using SharpMath;
using SharpMath.EquationsSystem.Preconditions;
using SharpMath.EquationsSystem.Solver;
using SharpMath.FiniteElement;
using SharpMath.FiniteElement._2D;
using SharpMath.FiniteElement._2D.Assembling;
using SharpMath.FiniteElement.Core.Assembling.Boundary.First;
using SharpMath.FiniteElement.Core.Assembling.Boundary.Second;
using SharpMath.FiniteElement.Core.Assembling.Params;
using SharpMath.FiniteElement.Core.Harmonic;
using SharpMath.FiniteElement.Core.Harmonic.Solution;
using SharpMath.FiniteElement.Materials.HarmonicWithoutChi;
using SharpMath.FiniteElement.Providers.Density;
using SharpMath.Geometry;
using SharpMath.Geometry._2D;
using SharpMath.Matrices.Converters;
using SharpMath.Matrices.Sparse;
using Vector = SharpMath.Vectors.Vector;

namespace InverseTask.DirectTask;

public class DirectSolver : IDirectSolver
{
    private readonly LocalOptimalSchemeConfig _losConfig;
    private HarmonicContext<Point, Element, SparseMatrix> _context = null!;
    private EquationAssembler _assembler = null!;

    private const double MagneticConstant = 4d * Math.PI * 1e-7;

    public DirectSolver(LocalOptimalSchemeConfig losConfig)
    {
        _losConfig = losConfig;
    }

    public void Allocate(Grid<Point, Element> grid)
    {
        _context = CreateContext(grid);
        _assembler = CreateAssembler(_context);
    }

    public double[] Solve(
        double frequency, 
        IMaterialProvider<Material> materialProvider, 
        Vector measuringPoints,
        double[] resultMemory
    )
    {
        _assembler.BuildEquation(_context)
            .ApplySecondConditions(_context)
        .ApplyFirstBoundary(_context);

        var profile = MatrixConverter.Convert(_context.Equation.Matrix);
        var equationProfile = new Equation<ProfileMatrix>(profile, _context.Equation.Solution, _context.Equation.RightSide);
        var profileSolver = new LUProfile();
        profileSolver.Solve(equationProfile);
        var solution = new ImpedanceSolution(_context.Grid, _context.Materials, equationProfile.Solution, frequency);

        foreach (var (x, pointIndex) in measuringPoints)
        {
            var impedance = solution.Calculate(new Point(x, 0));
            var roSeeming = impedance * impedance / (MagneticConstant * frequency);
            
            if (double.IsNaN(roSeeming))
            {
                throw new Exception("Вырожденное решение");
            }

            resultMemory[pointIndex] = roSeeming;
        }

        return resultMemory;
    }

    private HarmonicContext<Point, Element, SparseMatrix> CreateContext(Grid<Point, Element> grid)
    {
        var portraitBuilder = new MatrixPortraitBuilder();
        var matrix = portraitBuilder.Build(grid.Elements, grid.Nodes.Length);

        var context = new HarmonicContext<Point, Element, SparseMatrix>
        {
            Grid = grid,
            Equation = new Equation<SparseMatrix>(
                Matrix: matrix,
                RightSide: Vector.Create(grid.Nodes.Length * 2),
                Solution: Vector.Create(grid.Nodes.Length * 2)
            ),
            Materials = null,
            DensityFunctionProvider = null,
            FirstConditions = CreateFirst(grid),
            SecondConditions = CreateSecond(grid),
            Frequency = default
        };
        context.DensityFunctionProvider = new AnalyticComplexDensity(context, p => new Complex(0, 0));

        return context;
    }

    private FirstCondition[] CreateFirst(Grid<Point, Element> grid)
    {
        var result = new FirstCondition[grid.Nodes.XLength * 2];

        foreach (var element in EnumerateBottomElements(grid))
        {
            var nodeIndexes = element.GetBoundNodeIndexes(Bound.Bottom);

            foreach (var nodeIndex in nodeIndexes)
            {
                result[nodeIndex * 2] = new FirstCondition(nodeIndex * 2, 0);
                result[nodeIndex * 2 + 1] = new FirstCondition(nodeIndex * 2 + 1, 0);
            }
        }

        return result;
    }

    private SecondCondition[] CreateSecond(Grid<Point, Element> grid)
    {
        return EnumerateTopElementIndexes(grid)
            .Select(elementIndex => new SecondCondition(elementIndex, Bound.Bottom, [1d, 1d], ComponentType.Real))
            .ToArray();
    }

    private static IEnumerable<Element> EnumerateBottomElements(Grid<Point, Element> grid)
    {
        var xAxisLength = grid.Nodes.XLength;

        for (var i = 0; i < xAxisLength - 1; i++)
        {
            yield return grid.Elements[i];
        }
    }

    private static IEnumerable<int> EnumerateTopElementIndexes(Grid<Point, Element> grid)
    {
        var xAxisElements = grid.Nodes.XLength - 1;
        var yAxisElements = grid.Nodes.YLength - 1;

        for (var i = xAxisElements * (yAxisElements - 1); i <= xAxisElements * yAxisElements - 1; i++)
        {
            yield return i;
        }
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
}
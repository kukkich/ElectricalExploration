using SharpMath.FiniteElement._2D;
using SharpMath.FiniteElement._2D.Assembling;
using SharpMath.FiniteElement.Core.Assembling;
using SharpMath.Geometry;
using SharpMath.Geometry._2D;
using SharpMath.Matrices.Sparse;
using System;
using SharpMath;
using SharpMath.Vectors;
using SharpMath.FiniteElement.Materials.HarmonicWithoutChi;
using SharpMath.FiniteElement.Materials.Providers;
using SharpMath.FiniteElement;
using SharpMath.FiniteElement.Providers.Density;
using SharpMath.FiniteElement.Core.Assembling.Boundary.First;
using System.Reflection.Metadata.Ecma335;
using SharpMath.FiniteElement.Core.Assembling.Boundary.Second;
using SharpMath.FiniteElement.Core.Harmonic;
using SharpMath.EquationsSystem.Solver;
using SharpMath.FiniteElement.Core.Harmonic.Solution;
using SharpMath.Matrices.Converters;

namespace Harmonic2D.Tests;

public class ManualTest
{
    public void Run()
    {
        var grid = new Grid<Point, Element>(
            new PointsCollection(
                new [] {-1, 0d},
                new [] {-2 * 1e5, -1e5, -2 * 1e4, -1e4, 0}
            ),
            new List<Element>
            {
                new(new[] {0, 1, 5, 6}, 1e5, 1, 0),
                new(new[] {1, 2, 6, 7}, Math.Abs(-2e4 - (-1e5)), 1, 0),
                new(new[] {2, 3, 7, 8}, Math.Abs(-1e4 - (-2e4)), 1, 1),
                new(new[] {3, 4, 8, 9}, Math.Abs(0 - (-1e4)), 1, 0),
            }
        );

        var material = new Material[]
        {
            new (1d/(4*Math.PI * 1e-7), 0.01, 2*Math.PI*0.01),
            new (1d/(4*Math.PI * 1e-7), 0.1, 2*Math.PI*0.01),
        };

        IMatrixPortraitBuilder<SparseMatrix, Element> portraitBuilder = new MatrixPortraitBuilder();
        var matrix = portraitBuilder.Build(grid.Elements, grid.Nodes.Length);

        var equation = new Equation<SparseMatrix>(
            Matrix: matrix,
            RightSide: Vector.Create(grid.Nodes.Length * 2),
            Solution: Vector.Create(grid.Nodes.Length * 2)
        );

        var context = new Context<Point, Element, SparseMatrix>
        {
            Grid = grid,
            Equation = new Equation<SparseMatrix>(
                Matrix: matrix,
                RightSide: Vector.Create(grid.Nodes.Length * 2),
                Solution: Vector.Create(grid.Nodes.Length * 2)
            ),
            Materials = new FromArrayMaterialProvider(material),
            DensityFunctionProvider = null,
            FirstConditions = null,
            SecondConditions = null,
        };

        context.DensityFunctionProvider = new AnalyticComplexDensity(context, _ => 0);
        var first = new List<FirstCondition>()
        {
            new (0 * 2, 0),
            new (0 * 2 + 1, 0),
            new (5 * 2 + 1, 0),
            new (5 * 2 + 1, 0),
        };
        var second = new List<SecondCondition>()
        {
            new(3, Bound.Right, new[] {1d, 1d}, ComponentType.Real),
            new(3, Bound.Right, new[] {0d, 0d}, ComponentType.Imaginary)
        };

        context.FirstConditions = first.ToArray();
        context.SecondConditions = second.ToArray();

        var inserter = new Inserter();

        var assembler = new EquationAssembler(
            context,
            new LocalAssembler(context),
            inserter,
            new GaussExcluderSparse(),
            new SecondBoundaryApplier(context, inserter)
        );

        assembler.BuildEquation(context)
            .ApplySecondConditions(context)
            .ApplyFirstBoundary(context);

        var profile = MatrixConverter.Convert(context.Equation.Matrix);
        var equationProfile = new Equation<ProfileMatrix>(profile, context.Equation.Solution, context.Equation.RightSide);
        var profileSolver = new LUProfile();
        profileSolver.Solve(equationProfile);
        var solution = new FiniteElementSolution2DHarmonic(grid, context.Materials, equationProfile.Solution);
        Console.WriteLine(solution.Calculate(new Point(0, 0), 0));
        for (int i = 0; i < equationProfile.Solution.Length; i++)
        {
            equationProfile.Solution[i] *= 100;
        }
    }
}
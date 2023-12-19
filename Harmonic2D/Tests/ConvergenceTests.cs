using SharpMath.FiniteElement;
using SharpMath.FiniteElement._2D;
using SharpMath.FiniteElement.Core.Harmonic;
using SharpMath.FiniteElement.Materials.Providers;
using SharpMath.Geometry._2D;
using SharpMath.Geometry.Splitting;
using SharpMath.Matrices.Sparse;
using static System.Math;

namespace Harmonic2D.Tests;

public class ConvergenceTests
{
    public static Context<Point, Element, SparseMatrix> ExpProduct(int xSize, int ySize)
    {
        double UExpected(Point p, double t) => Exp(p.X * p.Y) * (Sin(t) + Cos(t));
        var testBuilder = new TestBuilder();
        var grid = new GridBuilder()
            .SetXAxis(new AxisSplitParameter(new double[]
            {
                0, 2
            }, new UniformSplitter(xSize)))
            .SetYAxis(new AxisSplitParameter(new double[]
            {
                0, 2
            }, new UniformSplitter(ySize)))
            .Build();

        var equation = testBuilder
            .SetFSin(p => -1d*Exp(p.X * p.X) * (1d + p.X * p.X + p.Y * p.Y))
            .SetFCos(p => Exp(p.X * p.X)*(1d - p.X*p.X - p.Y*p.Y))
            .SetUSin(p => Exp(p.X * p.X))
            .SetUCos(p => Exp(p.X * p.X))
            .SetMaterial(new DefaultMaterialProvider())
            .SetGrid(grid)
            .SetFirstAround(xSize, ySize, ComponentType.Imaginary)
            .SetFirstAround(xSize, ySize, ComponentType.Real)
            .SetSecondBoundary(Array.Empty<BoundaryCreation>())
            .Build();

        return equation;
    }
    public static Context<Point, Element, SparseMatrix> ExpSum(int xSize, int ySize)
    {
        double UExpected(Point p, double t) => Exp(p.X + p.Y) * Sin(t) + Exp(p.X - p.Y) * Cos(t);

        var testBuilder = new TestBuilder();
        var grid = new GridBuilder()
            .SetXAxis(new AxisSplitParameter(new double[]
            {
                0, 2
            }, new UniformSplitter(xSize)))
            .SetYAxis(new AxisSplitParameter(new double[]
            {
                0, 2
            }, new UniformSplitter(ySize)))
            .Build();

        var equation = testBuilder
            .SetFSin(p => -1d * Exp(p.X - p.Y) - 2d * Exp(p.X + p.Y))
            .SetFCos(p => Exp(p.X + p.Y) - 2d * Exp(p.X - p.Y))
            .SetUSin(p => Exp(p.X + p.Y))
            .SetUCos(p => Exp(p.X - p.Y))
            .SetMaterial(new DefaultMaterialProvider())
            .SetGrid(grid)
            .SetFirstAround(xSize, ySize, ComponentType.Imaginary)
            .SetFirstAround(xSize, ySize, ComponentType.Real)
            .SetSecondBoundary(Array.Empty<BoundaryCreation>())
            .Build();

        return equation;
    }
    public static Context<Point, Element, SparseMatrix> Linear()
    {
        double UExpected(Point p, double t) => (p.X + p.Y) * Sin(t) + (p.X - p.Y) * Cos(t);

        var testBuilder = new TestBuilder();
        var grid = new GridBuilder()
            .SetXAxis(new AxisSplitParameter(new double[]
            {
                0, 2
            }, new UniformSplitter(2)))
            .SetYAxis(new AxisSplitParameter(new double[]
            {
                0, 2
            }, new UniformSplitter(2)))
            .Build();

        var equation = testBuilder
            .SetFSin(p => 1d * (-1d * p.X + p.Y))
            .SetFCos(p => 1d * (p.X + p.Y))
            .SetUSin(p => p.X + p.Y)
            .SetUCos(p => p.X - p.Y)
            .SetMaterial(new DefaultMaterialProvider())
            .SetGrid(grid)
            .SetFirstAround(2, 2, ComponentType.Imaginary)
            .SetFirstAround(2, 2, ComponentType.Real)
            .SetSecondBoundary(Array.Empty<BoundaryCreation>())
            .Build();

        return equation;
    }

    public static Context<Point, Element, SparseMatrix> Square(int xSize, int ySize)
    {
        double UExpected(Point p, double t) => (p.X*p.X + p.Y*p.Y) * Sin(t) + (p.X * p.X - p.Y * p.Y) * Cos(t);

        var testBuilder = new TestBuilder();
        var grid = new GridBuilder()
            .SetXAxis(new AxisSplitParameter(new double[]
            {
                0, 2
            }, new UniformSplitter(xSize)))
            .SetYAxis(new AxisSplitParameter(new double[]
            {
                0, 2
            }, new UniformSplitter(ySize)))
            .Build();

        var equation = testBuilder
            .SetFSin(p => 1d * (-4d - p.X * p.X + p.Y * p.Y))
            .SetFCos(p => 1d * (p.X * p.X + p.Y * p.Y))
            .SetUSin(p => (p.X * p.X + p.Y * p.Y))
            .SetUCos(p => (p.X * p.X - p.Y * p.Y))
            .SetMaterial(new DefaultMaterialProvider())
            .SetGrid(grid)
            .SetFirstAround(xSize, ySize, ComponentType.Imaginary)
            .SetFirstAround(xSize, ySize, ComponentType.Real)
            .SetSecondBoundary(Array.Empty<BoundaryCreation>())
            .Build();

        return equation;
    }
    public static Context<Point, Element, SparseMatrix> Cubic(int xSize, int ySize)
    {
        double UExpected(Point p, double t) => (p.X * p.X * p.X + p.Y * p.Y * p.Y) * Sin(t) 
                                               + (p.X * p.X * p.X - p.Y * p.Y * p.Y) * Cos(t);

        var testBuilder = new TestBuilder();
        var grid = new GridBuilder()
            .SetXAxis(new AxisSplitParameter(new double[]
            {
                0, 2
            }, new UniformSplitter(xSize)))
            .SetYAxis(new AxisSplitParameter(new double[]
            {
                0, 2
            }, new UniformSplitter(ySize)))
            .Build();

        var equation = testBuilder
            .SetFSin(p => (-1d*p.X * p.X * p.X + p.Y * p.Y * p.Y) - 6 * (p.X + p.Y))
            .SetFCos(p => (p.X * p.X * p.X + p.Y * p.Y * p.Y) - 6*(p.X-p.Y))
            .SetUSin(p => (p.X * p.X * p.X + p.Y * p.Y * p.Y))
            .SetUCos(p => (p.X * p.X * p.X - p.Y * p.Y * p.Y))
            .SetMaterial(new DefaultMaterialProvider())
            .SetGrid(grid)
            .SetFirstAround(xSize, ySize, ComponentType.Imaginary)
            .SetFirstAround(xSize, ySize, ComponentType.Real)
            .SetSecondBoundary(Array.Empty<BoundaryCreation>())
            .Build();

        return equation;
    }
}
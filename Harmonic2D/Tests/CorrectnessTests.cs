using SharpMath.FiniteElement;
using SharpMath.FiniteElement._2D;
using SharpMath.FiniteElement.Core.Harmonic;
using SharpMath.FiniteElement.Materials.HarmonicWithoutChi;
using SharpMath.FiniteElement.Materials.Providers;
using SharpMath.Geometry._2D;
using SharpMath.Geometry.Splitting;
using SharpMath.Matrices.Sparse;

namespace Harmonic2D.Tests;

public static class CorrectnessTests
{
    public static Context<Point, Element, SparseMatrix> SecondBoundaryTest()
    {
        var testBuilder = new TestBuilder();
        var grid = new GridBuilder()
            .SetXAxis(new AxisSplitParameter(new double[]
            {
                0, 2
            }, new IIntervalSplitter[]
            {
                new UniformSplitter(1),
            }))
            .SetYAxis(new AxisSplitParameter(new double[]
            {
                1, 3
            }, new IIntervalSplitter[]
            {
                new UniformSplitter(1),
            }))
            .Build();

        var equation = testBuilder
            .SetFSin(p => -1d * p.X + p.Y)
            .SetFCos(p => p.X + p.Y)
            .SetUSin(p => p.X + p.Y)
            .SetUCos(p => p.X - p.Y)
            .SetMaterial(new DefaultMaterialProvider())
            .SetGrid(grid)
            .SetFirstBoundary(Array.Empty<BoundaryCreation>())
            .SetSecondBoundary(new BoundaryCreation[]
            {
                new(0, new[] {Bound.Bottom, Bound.Right, Bound.Top, Bound.Left}, ComponentType.Real),
                new(0, new[] {Bound.Bottom, Bound.Right, Bound.Top, Bound.Left}, ComponentType.Imaginary),
            })
            .SetUSinXDerivative(_ => 1)
            .SetUSinYDerivative(_ => 1)
            .SetUCosXDerivative(_ => 1)
            .SetUCosYDerivative(_ => -1d)
            .Build();

        return equation;
    }

    public static Context<Point, Element, SparseMatrix> FirstBoundaryTest()
    {
        var testBuilder = new TestBuilder();

        #region grid
        var grid = new GridBuilder()
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
            }))
            .Build();
        #endregion

        var equation = testBuilder
            .SetFSin(p => 1d*(-1d * p.X + p.Y))
            .SetFCos(p => 1d*(p.X + p.Y))
            .SetUSin(p => p.X + p.Y)
            .SetUCos(p => p.X - p.Y)
            .SetMaterial(new DefaultMaterialProvider())
            .SetGrid(grid)
            .SetFirstAround(2, 2, ComponentType.Imaginary)
            .SetFirstAround(2, 2, ComponentType.Real)
            .SetSecondBoundary(Array.Empty<BoundaryCreation>())
            .SetUSinXDerivative(_ => 1)
            .SetUSinYDerivative(_ => 1)
            .SetUCosXDerivative(_ => 1)
            .SetUCosYDerivative(_ => -1d)
            .Build();

        return equation;
    }

}
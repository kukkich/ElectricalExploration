using System.Numerics;
using SharpMath.Geometry;
using SharpMath.FiniteElement._2D;
using SharpMath.FiniteElement.Core.Assembling.Params;
using SharpMath.FiniteElement.Materials;
using SharpMath.FiniteElement.Materials.HarmonicWithoutChi;
using SharpMath.Geometry._2D;
using Vector = SharpMath.Vectors.Vector;

namespace InverseTask.DirectTask.Result;

public class ImpedanceSolution
{
    private readonly Grid<Point, Element> _grid;
    private readonly IMaterialProvider<Material> _materialProvider;
    private readonly Vector _weights;
    private readonly double _frequency;

    public ImpedanceSolution(
        Grid<Point, Element> grid,
        IMaterialProvider<Material> materialProvider,
        Vector weights,
        double frequency
    )
    {
        _grid = grid;
        _materialProvider = materialProvider;
        _weights = weights;
        _frequency = frequency;
    }

    public double Calculate(Point point)
    {
        var element = _grid.Elements.First(x => ElementHas(x, point));
        var material = _materialProvider.GetById(element.MaterialId);

        var u = GetFunctionComponents(point, element);
        var electricity = _frequency * u.Magnitude;

        var uDerivative = GetDerivativeFunctionComponents(point, element);
        var magnetic = material.Lambda * uDerivative.Magnitude;
        var result = electricity / magnetic;

        return result;
    }

    private Complex GetFunctionComponents(Point point, IFiniteElement element)
    {
        var leftBottom = _grid.Nodes[element.NodeIndexes[0]];
        var rightTop = _grid.Nodes[element.NodeIndexes[^1]];
        var (xMin, xMax) = (leftBottom.X, rightTop.X);
        var (yMin, yMax) = (leftBottom.Y, rightTop.Y);
        var (hx, hy) = (xMax - xMin, yMax - yMin);

        Span<double> x = stackalloc double[2];
        Span<double> y = stackalloc double[2];

        x[0] = (xMax - point.X) / hx;
        y[0] = (yMax - point.Y) / hy;
        x[1] = (point.X - xMin) / hx;
        y[1] = (point.Y - yMin) / hy;

        Span<double> basicValues = stackalloc double[x.Length * y.Length];

        basicValues[0] = x[0] * y[0];
        basicValues[1] = x[1] * y[0];
        basicValues[2] = x[0] * y[1];
        basicValues[3] = x[1] * y[1];

        var (us, uc) = (0d, 0d);

        for (var i = 0; i < 4; i++)
        {
            var nodeIndex = element.NodeIndexes[i];
            us += _weights[nodeIndex * 2] * basicValues[i];
            uc += _weights[nodeIndex * 2 + 1] * basicValues[i];
        }

        return new Complex(us, uc);
    }

    private Complex GetDerivativeFunctionComponents(Point point, IFiniteElement element)
    {
        var leftBottom = _grid.Nodes[element.NodeIndexes[0]];
        var rightTop = _grid.Nodes[element.NodeIndexes[^1]];
        var (xMin, xMax) = (leftBottom.X, rightTop.X);
        var (yMin, yMax) = (leftBottom.Y, rightTop.Y);
        var (hx, hy) = (xMax - xMin, yMax - yMin);

        Span<double> x = stackalloc double[2];
        Span<double> y = stackalloc double[2];

        x[0] = (xMax - point.X) / hx;
        y[0] = -1d / hy;
        x[1] = (point.X - xMin) / hx;
        y[1] = 1d / hy;

        Span<double> basicValues = stackalloc double[x.Length * y.Length];

        basicValues[0] = x[0] * y[0];
        basicValues[1] = x[1] * y[0];
        basicValues[2] = x[0] * y[1];
        basicValues[3] = x[1] * y[1];

        var (us, uc) = (0d, 0d);

        for (var i = 0; i < 4; i++)
        {
            var nodeIndex = element.NodeIndexes[i];
            us += _weights[nodeIndex * 2] * basicValues[i];
            uc += _weights[nodeIndex * 2 + 1] * basicValues[i];
        }

        return new Complex(us, uc);
    }

    private bool ElementHas(IFiniteElement element, Point node)
    {
        var leftBottom = _grid.Nodes[element.NodeIndexes[0]];
        var rightTop = _grid.Nodes[element.NodeIndexes[^1]];

        return leftBottom.X <= node.X && node.X <= rightTop.X &&
               leftBottom.Y <= node.Y && node.Y <= rightTop.Y;
    }
}
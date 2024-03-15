using System.Collections;
using SharpMath;
using SharpMath.Matrices.Sparse.Storages;
using SharpMath.Vectors;

namespace InverseTask.EquationSystem;

public class ComputedVector : IReadonlyVector<double>
{
    public double this[int x] => _f[x] - _alpha[x] * (_parameterCurrent[x] - _parameterInitial[x]);

    public int Length => _alpha.Length;
    public double Norm => Math.Sqrt(Vector.ScalarProduct(this, this));

    private readonly Vector _parameterCurrent;
    private readonly Vector _parameterInitial;
    private readonly Vector _alpha;
    private readonly Vector _f;

    public ComputedVector(ParameterDirectionEquation equation)
    {
        _parameterCurrent = equation.ParameterCurrent;
        _parameterInitial = equation.ParameterInitial;
        _alpha = equation.Alpha;
        _f = equation.F;
    }

    public double ScalarProduct(IReadonlyVector<double> v)
    {
        return Vector.ScalarProduct(this, v);
    }

    public IEnumerator<IndexValue<double>> GetEnumerator()
    {
        for (var i = 0; i < Length; i++)
            yield return new IndexValue<double>(this[i], i);
    }

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }
}
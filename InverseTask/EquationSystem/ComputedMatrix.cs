using SharpMath.Matrices;
using SharpMath.Vectors;

namespace InverseTask.EquationSystem;

public class ComputedMatrix : MatrixBase
{
    public override double this[int row, int column]
    {
        get
        {
            if (row != column)
            {
                return _matrix[row, column];
            }

            return _matrix[row, column] + _alpha[row];
        }
    }

    public override int Rows => _matrix.Rows;
    public override int Columns => _matrix.Columns;

    private readonly Matrix _matrix;
    private readonly Vector _alpha;

    public ComputedMatrix(ParameterDirectionEquation equation)
    {
        _matrix = equation.A;
        _alpha = equation.Alpha;
    }

    public override ImmutableMatrix AsImmutable()
    {
        throw new NotImplementedException();
    }

    public override Matrix AsMutable()
    {
        throw new NotImplementedException();
    }
}
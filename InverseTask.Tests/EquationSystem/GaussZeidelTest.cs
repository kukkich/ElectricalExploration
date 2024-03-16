using InverseTask.EquationSystem;
using Microsoft.Extensions.Logging.Abstractions;
using SharpMath;
using SharpMath.Matrices;
using SharpMath.Matrices.Transformation;
using SharpMath.Vectors;

namespace InverseTask.Tests.EquationSystem;

public class GaussZeidelTest
{
    private GaussZeidelSolver Solver = null!;
    private GaussZeidelConfig Config = null!;
    
    [SetUp]
    public void Setup()
    {
        Config = new GaussZeidelConfig()
        {
            Precision = 1e-8,
            MaxIteration = 100,
            Relaxation = 1
        };

        Solver = new GaussZeidelSolver(Config, NullLogger.Instance);
    }

    [Test]
    public void GIVEN_slae_WHEN_solve_THEN_correct_precision_expected()
    {
        var matrix = new Matrix(new double[,]
        {
            { 5, -1, 2 },
            { 1, 7, -3 },
            { -3, 2, 5 },
        });
        var rightSide = new Vector(9, 6, 16);

        Solver.Allocate(rightSide.Length);
        var actualSolution = Solver.Solve(matrix, rightSide);

        var relativeDiscrepancy = LinAl.Subtract(
                                      rightSide,
                                      LinAl.Multiply(
                                          matrix, 
                                          actualSolution
                                      )
                                  ).Norm / rightSide.Norm;
        
        Assert.That(relativeDiscrepancy, Is.LessThanOrEqualTo(Config.Precision));
    }
}
import org.example.Gen;
import org.example.HouseholderInverse;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;

public class GenTest {

    @Test
    public void GenTest1() {
        int n = 100;
        double alpha = 1;
        double beta = 2;
        double[][] a = new double[n][n];
        double[][] a_inv = new double[n][n];

        Gen g = new Gen();
        g.mygen(a, a_inv, n, alpha, beta, 1, 2, 2, 1);

        double[][] a_copy = new double[n][n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                a_copy[i][j] = a[i][j];

        HouseholderInverse.invert(a_copy);

        System.out.println("Original A:");
        g.print_matr(a, n);
        System.out.println("Exact inverse (generator):");
        g.print_matr(a_inv, n);
        System.out.println("Your inverse:");
        g.print_matr(a_copy, n);

        double[][] r = new double[n][n];
        g.matr_mul(a, a_copy, r, n);
        for (int i = 0; i < n; i++)
            r[i][i] -= 1.0;

        double residual = g.matr_inf_norm(r, n);
        System.out.println("Residual norm: " + residual);

        double TOLERANCE = 1e-10;
        assertTrue(residual < TOLERANCE);
    }
}
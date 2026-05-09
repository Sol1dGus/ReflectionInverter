import org.example.Gen;
import org.example.ReflectionInverse;
import org.junit.jupiter.api.Test;

public class GenTest {

    /**
     * Тест 1: Исследуем рост погрешности при увеличении beta (beta >> alpha).
     * alpha фиксирован = 1, beta меняется логарифмически от 10^1 до 10^8.
     */
    @Test
    public void generateMatrixTable1() {
        System.out.println("Таблица 1: Матрица простой структуры. beta > alpha (alpha = 1, beta от 10^1 до 10^8)\n");
        generateMatrixTable(100, 1.0, 1.0, 1.0, 10.0, 8);
    }

    /**
     * Тест 2: Исследуем рост погрешности при уменьшении alpha (alpha << beta).
     * beta фиксирован = 1, alpha меняется логарифмически от 10^-1 до 10^-8.
     */
    @Test
    public void generateMatrixTable2() {
        System.out.println("Таблица 2: Матрица простой структуры. alpha < 1 (alpha от 10^-1 до 10^-8, beta = 1)\n");
        generateMatrixTable(100, 1.0,  1.0, 0.1, 1.0, 8);
    }

    private void generateMatrixTable(int n, double alpha_min, double beta_min,  double a_coef, double b_coef, int num_points) {
        // Печать заголовка таблицы
        System.out.println("+------------+------------+--------------+----------------+--------------+--------------+--------------+--------------+");
        System.out.println("|     alpha     |    beta       |     ||A||inf      |  ||A^{-1}||inf    |     cond(A)    |    ||z||inf       |     zeta         |      ||r||inf     |");
        System.out.println("+------------+------------+--------------+----------------+--------------+--------------+--------------+--------------+");
        double alpha, beta;
        alpha = alpha_min;
        beta = beta_min;

        for (int i = 0; i <= num_points; i++) {
            // Коэффициент интерполяции [0, 1]
            double t = (double) Math.max(0, i) / num_points;

            // Вычисление текущих параметров
            alpha = alpha * a_coef;
            beta = beta * b_coef;

            try {
                // Анализ для текущей пары (alpha, beta)
                analyzeSimpleStructure(n, alpha, beta, 0, 1);
            } catch (Exception e) {
                // Обработка ошибок
                System.out.printf("|%12.4e|%12.4e|   ERROR      |   ERROR        |   ERROR      |   ERROR      |   ERROR      |   ERROR      |   ERROR      |\n", alpha, beta);
            }
        }
        System.out.println("+------------+------------+--------------+----------------+--------------+--------------+--------------+--------------+\n");
    }

    /**
     * Основной цикл анализа: генерация, обращение, валидация через тестовую СЛАУ.
     */
    private void analyzeSimpleStructure(int n, double alpha, double beta, int variant, int schema) {
        Gen wrapper = new Gen();
        double[][] a = new double[n][n];           // Исходная матрица
        double[][] a_inv_exact = new double[n][n]; // Точная обратная

        // Генерация тестовой матрицы с заданными параметрами
        wrapper.mygen(a, a_inv_exact, n, alpha, beta, 0, 2, variant, schema);

        // Копия для обращения (исходная a сохраняется для вычисления невязки)
        double[][] a_copy = new double[n][n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                a_copy[i][j] = a[i][j];

        // Вычисление обратной матрицы методом отражений (in-place)
        ReflectionInverse.invert(a_copy);

        double[][] z = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                z[i][j] = a_copy[i][j] - a_inv_exact[i][j];
            }
        }

        // Норма матрицы ошибки
        double norm_z = wrapper.matr_inf_norm(z, n);

        // Матрица невязки R
        double[][] r = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double sum = 0.0;
                for (int k = 0; k < n; k++) {
                    sum += a[i][k] * a_copy[k][j];
                }
                r[i][j] = (i == j) ? sum - 1.0 : sum;
            }
        }

        // Норма матрицы невязки
        double norm_r = wrapper.matr_inf_norm(r, n);

        // 6. Нормы матриц и число обусловленности
        double norm_A = wrapper.matr_inf_norm(a, n);
        double norm_Ainv_star = wrapper.matr_inf_norm(a_inv_exact, n);
        double norm_Ainv = wrapper.matr_inf_norm(a_copy, n); // Норма вычисленной обратной
        double cond = norm_A * norm_Ainv;

        // 7. Относительные метрики
        double rel_error = norm_z / norm_Ainv_star;

        // Вывод строки результатов в таблицу
        System.out.printf("|%12.2e|%12.2e|%14.2e|%16.2e|%14.2e|%14.4e|%14.4e|%14.4e|\n",
                alpha, beta, norm_A, norm_Ainv, cond, norm_z, rel_error, norm_r);
    }
}
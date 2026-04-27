import org.example.Gen;
import org.example.ReflectionInverse;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;

public class GenTest {

    /**
     * Тест 1: Исследуем рост погрешности при увеличении beta (beta >> alpha).
     * alpha фиксирован = 1, beta меняется логарифмически от 10^1 до 10^8.
     * Ожидаем: рост cond(A) → рост ошибки решения.
     */
    @Test
    public void generateMatrixTable1() {
        System.out.println("Таблица 1: Матрица простой структуры. beta > alpha (alpha = 1, beta от 10^1 до 10^8)\n");
        generateMatrixTable(100, 1.0, 1.0, 10.0, 1e8, 8, true);
    }

    /**
     * Тест 2: Исследуем рост погрешности при уменьшении alpha (alpha << beta).
     * beta фиксирован = 1, alpha меняется логарифмически от 10^-1 до 10^-8.
     * Ожидаем: аналогичный рост обусловленности и ошибки.
     */
    @Test
    public void generateMatrixTable2() {
        System.out.println("Таблица 2: Матрица простой структуры. alpha < 1 (alpha от 10^-1 до 10^-8, beta = 1)\n");
        generateMatrixTable(100, 1e-1, 1e-8, 1.0, 1.0, 7, true);
    }

    /**
     * Генерирует сетку параметров (alpha, beta) и выводит сводную таблицу метрик.
     * @param n размер матрицы
     * @param alpha_min/max границы изменения alpha
     * @param beta_min/max границы изменения beta
     * @param num_points количество точек в сетке
     * @param log_scale если true — логарифмическая сетка, иначе линейная
     */
    private void generateMatrixTable(int n, double alpha_min, double alpha_max,
                                     double beta_min, double beta_max,
                                     int num_points, boolean log_scale) {
        // Печать заголовка таблицы
        System.out.println("+------------+------------+--------------+----------------+--------------+--------------+--------------+--------------+--------------+");
        System.out.println("|    alpha   |    beta    |   ||A||inf   |  ||A^{-1}||inf |   cond(A)    |    ||z||inf  |     zeta     |    ||r||inf  |     rho      |");
        System.out.println("+------------+------------+--------------+----------------+--------------+--------------+--------------+--------------+--------------+");

        for (int i = 0; i <= num_points; i++) {
            // Коэффициент интерполяции [0, 1]
            double t = (double) Math.max(0, i) / num_points;
            double alpha, beta;

            // Вычисление текущих параметров
            if (log_scale) {
                alpha = alpha_min * Math.pow(alpha_max / alpha_min, t);
                beta = beta_min * Math.pow(beta_max / beta_min, t);
            } else {
                alpha = alpha_min + t * (alpha_max - alpha_min);
                beta = beta_min + t * (beta_max - beta_min);
            }

            try {
                // Анализ для текущей пары (alpha, beta)
                analyzeSimpleStructure(n, alpha, beta, 1, 1);
            } catch (Exception e) {
                // Обработка ошибок
                System.out.printf("|%12.4e|%12.4e|   ERROR      |   ERROR        |   ERROR      |   ERROR      |   ERROR      |   ERROR      |   ERROR      |\n", alpha, beta);
            }
        }
        System.out.println("+------------+------------+--------------+----------------+--------------+--------------+--------------+--------------+--------------+\n");
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

        // --- Валидация через решение тестовой СЛАУ ---

        // 1. Формируем эталонный вектор решения x*
        double[] x_star = new double[n];
        for (int i = 0; i < n; i++) {
            x_star[i] = Math.sin(i + 1);
        }

        // 2. Вычисляем правую часть: f = A * x*
        double[] f = new double[n];
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                sum += a[i][j] * x_star[j];
            }
            f[i] = sum;
        }

        // 3. Решаем систему через ВЫЧИСЛЕННУЮ обратную: x~ = A^{-1}_comp * f
        double[] x_tilde = new double[n];
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                sum += a_copy[i][j] * f[j]; // a_copy теперь содержит A^{-1}
            }
            x_tilde[i] = sum;
        }

        // 4. Вычисляем абсолютную ошибку решения: z = x~ - x*
        double norm_z = 0.0;
        double norm_x_star = 0.0;
        for (int i = 0; i < n; i++) {
            norm_z = Math.max(norm_z, Math.abs(x_tilde[i] - x_star[i]));
            norm_x_star = Math.max(norm_x_star, Math.abs(x_star[i]));
        }

        // 5. Вычисляем невязку: r = A*x~ - f
        double norm_r = 0.0;
        double norm_f = 0.0;
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                sum += a[i][j] * x_tilde[j]; // Используем исходную A
            }
            norm_r = Math.max(norm_r, Math.abs(sum - f[i]));
            norm_f = Math.max(norm_f, Math.abs(f[i]));
        }

        // 6. Нормы матриц и число обусловленности
        double norm_A = wrapper.matr_inf_norm(a, n);
        double norm_Ainv = wrapper.matr_inf_norm(a_copy, n); // Норма вычисленной обратной
        double cond = norm_A * norm_Ainv;

        // 7. Относительные метрики
        double rel_error = norm_z / norm_x_star;   // zeta: относительная ошибка решения
        double rel_residual = norm_r / norm_f;     // rho: относительная невязка

        // Вывод строки результатов в таблицу
        System.out.printf("|%12.2e|%12.2e|%14.2e|%16.2e|%14.2e|%14.4e|%14.4e|%14.4e|%14.4e|\n",
                alpha, beta, norm_A, norm_Ainv, cond, norm_z, rel_error, norm_r, rel_residual);
    }
}
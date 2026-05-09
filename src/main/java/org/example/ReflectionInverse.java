package org.example;

/**
 * Обращение матрицы методом отражений.
 * Схема: A = Q·R → A⁻¹ = R⁻¹·Qᵀ
 */
public class ReflectionInverse {
    private static final double EPS = 1e-14;

    /**
     * Точка входа: обращает матрицу A на месте.
     * После выполнения a содержит A⁻¹.
     */
    public static void invert(double[][] a) {
        int n = a.length;
        if (n == 0) return;

        qrDecompositionInPlace(a);
        invertUpperTriangularInPlace(a);
        multiplyByQTInPlace(a);
    }

    // =========================
    // ЭТАП 1: QR-разложение
    // =========================
    /**
     * Приводит A к верхнетреугольному виду R методом отражений.
     * In-place схема:
     *   - верхний треугольник → элементы R
     *   - строгий нижний треугольник → векторы отражений w
     */
    private static void qrDecompositionInPlace(double[][] a) {
        int n = a.length;

        for (int j = 0; j < n - 1; j++) {
            double norm = colNorm(a, j, j);
            if (norm < EPS)
                throw new RuntimeException("Matrix is singular");

            // Выбор знака для численной устойчивости
            double alpha = (a[j][j] >= 0) ? -norm : norm;
            double v0 = a[j][j] - alpha;

            if (Math.abs(v0) < EPS) {
                a[j][j] = alpha;
                continue;
            }

            for (int i = j + 1; i < n; i++)
                a[i][j] /= v0;

            a[j][j] = alpha;

            double sumsq = colSumsq(a, j, j + 1);
            double wTw = 1.0 + sumsq;
            double beta = 2.0 / wTw;

            for (int k = j + 1; k < n; k++) {
                double dot = a[j][k];
                for (int i = j + 1; i < n; i++)
                    dot += a[i][j] * a[i][k];

                a[j][k] -= beta * dot;
                for (int i = j + 1; i < n; i++)
                    a[i][k] -= beta * dot * a[i][j];
            }
        }
    }

    // =========================
    // ЭТАП 2: Инверсия верхней треугольной матрицы R
    // =========================
    /**
     * Заменяет верхний треугольник R на R⁻¹ (обратная тоже верхнетреугольная).
     * Формула: R⁻¹[i,i] = 1/R[i,i]; R⁻¹[j,i] = -Σ(k=j+1..i) R[j,k]·R⁻¹[k,i] / R[j,j]
     */
    private static void invertUpperTriangularInPlace(double[][] a) {
        int n = a.length;

        for (int i = n - 1; i >= 0; i--) {
            if (Math.abs(a[i][i]) < EPS)
                throw new RuntimeException("Matrix is singular");

            double rii = a[i][i];
            a[i][i] = 1.0 / rii;

            for (int j = i - 1; j >= 0; j--) {
                double sum = 0.0;
                for (int k = j + 1; k <= i; k++)
                    sum += a[j][k] * a[k][i];

                a[j][i] = -sum / a[j][j];
            }
        }
    }

    // =========================
    // ЭТАП 3: Умножение R⁻¹ на Qᵀ справа
    // =========================
    /**
     * Применяет накопленные отражения в обратном порядке: результат = R⁻¹ · H₁·H₂·...·Hₙ₋₁ = R⁻¹·Qᵀ.
     * Векторы отражений извлекаются из нижнего треугольника матрицы.
     */
    private static void multiplyByQTInPlace(double[][] a) {
        int n = a.length;
        for (int j = n - 2; j >= 0; j--) {
            int len = n - j; // Размер вектора

            double[] w = new double[len];
            w[0] = 1.0;
            for (int i = 1; i < len; i++)
                w[i] = a[j + i][j];

            double wTw = vecSumsq(w);
            double beta = 2.0 / wTw;

            // Очистка нижнего треугольника
            for (int i = j + 1; i < n; i++)
                a[i][j] = 0.0;

            double[] rowCopy = new double[len]; // Буфер для строки

            for (int i = 0; i < n; i++) {
                System.arraycopy(a[i], j, rowCopy, 0, len);

                double dot = 0.0;
                for (int k = 0; k < len; k++)
                    dot += rowCopy[k] * w[k];

                for (int k = 0; k < len; k++)
                    a[i][j + k] = rowCopy[k] - beta * dot * w[k];
            }
        }
    }

    // ---------------------- вспомогательные методы ----------------------

    /** Евклидова норма хвоста столбца  */
    private static double colNorm(double[][] a, int col, int start) {
        int n = a.length;
        double r = 0.0;
        for (int i = start; i < n; i++) r = Math.hypot(r, a[i][col]);
        return r;
    }

    /** Сумма квадратов элементов хвоста столбца */
    private static double colSumsq(double[][] a, int col, int start) {
        int n = a.length;
        double s = 0.0;
        for (int i = start; i < n; i++) s += a[i][col] * a[i][col];
        return s;
    }

    /** Сумма квадратов элементов вектора */
    private static double vecSumsq(double[] v) {
        double s = 0.0;
        for (double x : v) s += x * x;
        return s;
    }
}
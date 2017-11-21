package venherak.com.lab5;

import static java.lang.Math.abs;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import static java.lang.System.arraycopy;

public class Model {
    private static String txt = "";
    private static String planMatrix = "";
    private static String regresEquation = "";
    private static String student = "";
    private static String fisher = "";
    private static String naturalizetEquation = "";
    private static String findings = "";

    private double[][] y;
    private int numberOfCoeff = 4;
    private int k = 3; //Колонки
    private int N = 4; //Строки
    private int p = 0; //Дробность
    private double l = 0; //звёздное плечо

    private double[] yser;
    private double[] disp;

    private double Gp = 0;
    private double GTable;
    private double STable;

    private int[] f;
    private double Fp;
    private double[] Akoef;
    private double sb;
    private double[] b;
    private double[] t;
    private double[] fish;
    private double[] tmod;
    private double[] xDash2;
    private double[] xMin = {-6, 0, -6};
    private double[] xMax = {2, 2, 8};
    private double[] x0;
    private double[] xDel;
    private double yimax;
    private double yimin;

    private double[][] xMatrix;
    private double my;
    private static int m = 2;

    public Model(double x1min, double x1max, double x2min, double x2max, double x3min, double x3max) {
        xMin = new double[3];
        xMax = new double[3];

        xMin[0] = x1min;
        xMin[1] = x2min;
        xMin[2] = x3min;
        xMax[0] = x1max;
        xMax[1] = x2max;
        xMax[2] = x3max;
        x0 = new double[3];
        xDel = new double[3];
        x0[0] = (xMax[0] + xMin[0]) / 2;
        x0[1] = (xMax[1] + xMin[1]) / 2;
        x0[2] = (xMax[2] + xMin[2]) / 2;

        xDel[0] = abs(xMax[0] - xMin[0]) / 2;
        xDel[1] = abs(xMax[0] - xMin[0]) / 2;
        xDel[2] = abs(xMax[0] - xMin[0]) / 2;
    }

    public void calc() {
        int step = 1;
        do {
            step++;
            changeConditionsOfExperiment(step);
            m = 1;
            do {
                m++;
                determineResponseFunc();
                determineRegressionCoeff();
            } while ((!determineCohren()) && (m < 5));

            determineStudent();
        } while ((!determineFisher()) && (step < 2));
    }

    private void determineResponseFunc() {
        l = getL(p);
        refreshMatrix();

        xMatrix = new double[N][numberOfCoeff - 1];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < k; j++)
                xMatrix[i][j] = x0[j] + plan_matrix[i][j] * xDel[j];

            for (int j = 0; j < k - 1; j++)
                for (int q = j + 1; q < k; q++)
                    xMatrix[i][j + q + k - 1] = xMatrix[i][j] * xMatrix[i][q];
            xMatrix[i][6] = xMatrix[i][0] * xMatrix[i][1] * xMatrix[i][2];
            for (int j = 7; j < 10; j++)
                xMatrix[i][j] = pow(xMatrix[i][j - 7], 2);
        }
        f = new int[4];
        f[0] = m - 1;
        f[1] = N;
        f[2] = N * (m - 1);
        xDash2 = new double[k];
        for (int i = 0; i < k; i++) {
            xDash2[i] = (pow(2, k) + 2 * pow(l, 2)) / N;
        }
        for (int i = 0; i < N - 1; i++)
            for (int j = 7; j < 10; j++) {
                plan_matrix[i][j] -= xDash2[j - 7];
            }
        y = new double[N][m];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < m; j++) {
                y[i][j] = f(xMatrix[i][0], xMatrix[i][1], xMatrix[i][2])
                        + 10 * Math.random() - 5;
            }
        }
    }

    private void determineRegressionCoeff() {
        yser = new double[N];
        b = new double[numberOfCoeff];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < m; j++) {
                yser[i] += y[i][j];
            }
            yser[i] /= m;
        }
        my = 0;
        for (double aYser : yser) {
            my += aYser;
        }
        my /= yser.length;

        if (N <= 8) {
            double[][] a = new double[N - 1][N - 1];
            double[] A = new double[N - 1];

            double[] mx = new double[N - 1];

            for (int i = 0; i < N - 1; i++) {
                for (int j = 0; j < numberOfCoeff; j++) {
                    mx[i] += xMatrix[j][i];
                    A[i] += xMatrix[j][i] * yser[j];

                }
                mx[i] /= numberOfCoeff;
                A[i] /= numberOfCoeff;

                for (int q = 0; q < numberOfCoeff - 1; q++) {
                    for (int j = 0; j < N; j++)
                        a[i][q] += xMatrix[j][i] * xMatrix[j][q];
                    a[i][q] /= N;
                }
            }
            double[][] t1 = new double[numberOfCoeff][numberOfCoeff];
            double[][] t2 = new double[numberOfCoeff][numberOfCoeff];
            double[][] columns = new double[numberOfCoeff][numberOfCoeff];
            for (int i = 0; i < numberOfCoeff - 1; i++) {
                t1[i + 1][0] = mx[i];
                t2[i + 1][0] = mx[i];
                columns[0][i + 1] = mx[i];
                t1[0][i + 1] = mx[i];
                t2[0][i + 1] = mx[i];
                columns[i + 1][0] = mx[i];
                for (int j = 0; j < numberOfCoeff - 1; j++) {
                    t1[i + 1][j + 1] = a[i][j];
                    t2[i + 1][j + 1] = a[i][j];
                    columns[j + 1][i + 1] = a[i][j];
                }
            }
            t1[0][0] = 1;
            t2[0][0] = 1;
            columns[0][0] = 1;
            double t2det = determinantOfMatrix(t2);
            double[] row1 = new double[N];
            row1[0] = my;
            arraycopy(A, 0, row1, 1, N - 1);
            replaceColumn(0, t1, row1);
            b[0] = determinantOfMatrix(t1) / t2det;
            for (int i = 1; i < numberOfCoeff; i++) {
                replaceColumn(i - 1, t1, columns[i - 1]);
                replaceColumn(i, t1, row1);
                b[i] = determinantOfMatrix(t1) / t2det;
            }
        } else {
            for (int i = 0; i < numberOfCoeff; i++)
                b[i] = 0;
            b[0] = my;

            for (int i = 0; i < N; i++) {
                for (int j = 0; j < k; j++)
                    b[j + 1] += plan_matrix[i][j] * yser[i];
                for (int j = 0; j < k - 1; j++)
                    for (int q = j + 1; q < k; q++)
                        b[j + q + k] += plan_matrix[i][k - 1 + q + j] * yser[i];
                b[7] += plan_matrix[i][6] * yser[i];
                for (int j = 0; j < k; j++)
                    b[(int) pow(2, k) + j] += plan_matrix[i][(int) pow(2, k) + j - 1] * yser[i];
            }
            for (int i = 1; i < k + 1; i++) {
                b[i] /= pow(2, k - p) + 2 * l * l;
            }
            for (int i = k + 1; i < pow(2, k); i++) {
                b[i] /= (int) pow(2, k - p);
            }
            for (int i = (int) pow(2, k); i < numberOfCoeff; i++) {
                b[i] /= 2 * pow(l, 4);
                b[i] = sqrt(b[i]);
            }
            for (int i = 0; i < k; i++)
                b[0] -= xDash2[i] * b[(int) pow(2, k) + i];
        }
    }

    private boolean determineCohren() {
        disp = new double[N];
        double[] ts = new double[m];
        for (int i = 0; i < N; i++) {
            arraycopy(y[i], 0, ts, 0, m);
            disp[i] = getDispersion(ts);
        }
        double max = 0;
        double summ = 0;
        for (double t : disp) {
            if (max < t) {
                max = t;
            }
            summ += t;
        }
        Gp = max / summ;
        GTable = getCochrenTableValue(f[0], f[1]);
        return Gp < GTable;
    }

    private void determineStudent() {
        t = new double[numberOfCoeff];
        tmod = new double[numberOfCoeff];
        sb = 0;
        double[] beta = new double[numberOfCoeff];
        for (int i = 0; i < N; i++)
            sb += disp[i];
        sb /= N;
        beta[0] = my;
        for (int i = 1; i < numberOfCoeff; i++) {
            for (int j = 0; j < N; j++)
                beta[i] += xMatrix[j][i - 1] * yser[j];
            beta[i] /= N;
        }
        int st3 = (m - 1) * N;
        if (st3 > 30)
            st3 = 30;
        double sbeta = sqrt(sb / (N * m));
        if (N <= 8) {
            for (int i = 0; i < numberOfCoeff; i++) {
                tmod[i] = abs(beta[i]) / sbeta;
                if (tmod[i] > tableStudent[st3])
                    t[i] = 1;
                else
                    t[i] = 0;
            }
        } else {
            double[] sB = new double[numberOfCoeff];
            sB[0] = sb / (m * (N + 2 * pow(l, 4)));
            for (int i = 1; i < k + 1; i++) {
                sB[i] = sb / (m * (pow(2, k - p) + 2 * l * l));
            }
            for (int i = k + 1; i < pow(2, k); i++) {
                sB[i] = sb / (m * pow(2, k - p));
            }
            for (int i = (int) pow(2, k); i < numberOfCoeff; i++) {
                sB[i] = sb / (m * 2 * pow(l, 4));
            }
            STable = tableStudent[st3];
            for (int i = 0; i < numberOfCoeff; i++) {
                tmod[i] = abs(beta[i] / sbeta);
                if (tmod[i] > STable)
                    t[i] = 1;
                else
                    t[i] = 0;
            }
        }
    }

    private boolean determineFisher() {
        fish = new double[numberOfCoeff];
        for (int i = 0; i < numberOfCoeff; i++) {
            fish[i] += b[0] * t[0];
            for (int j = 1; j < numberOfCoeff; j++)
                if (!((N > 8) && ((j == 4) || (j == 6))))
                    fish[i] += b[j] * t[j] * plan_matrix[i][j - 1];
        }
        int d = 0;
        for (double aT : t)
            if (aT == 1)
                ++d;
        f[3] = N - d;
        Fp = 0;
        double Sad = 0;
        if (f[3] > 0) {
            for (int i = 0; i < numberOfCoeff; i++)
                Sad += pow((fish[i] - yser[i]), 2);
            Sad *= m;
            Sad /= f[3];
            Fp = Sad / sb;
        }
        return Fp < getFisherTableValue(f[2], f[3]);
    }

    private void normalizeBCoeff() {
        Akoef = new double[numberOfCoeff];
        Akoef[0] = b[0];
        for (int i = 0; i < k; i++) {
            Akoef[0] += sqrt(b[8 + i]) * pow(x0[i] / xDel[i], 2)
                    - b[i + 1] * x0[i] / xDel[i];
        }
        Akoef[0] += b[4] * x0[0] * x0[1] / (xDel[0] * xDel[1]) + b[6] * x0[1] * x0[2] / (xDel[1] * xDel[2]);
        Akoef[1] = b[1] / xDel[0] - 2 * b[8] * x0[0] / pow(xDel[0], 2) - b[4] * x0[1] / (xDel[0] * xDel[1]);
        Akoef[2] = b[2] / xDel[1] - 2 * b[9] * x0[1] / pow(xDel[1], 2) - b[4] * x0[0] / (xDel[0] * xDel[1]) - b[6] * x0[0] / (xDel[1] * xDel[2]);
        Akoef[3] = b[3] / xDel[2] - 2 * b[10] * x0[2] / pow(xDel[2], 2) - b[6] * x0[1] / (xDel[1] * xDel[2]);
        Akoef[4] = b[4] / (xDel[0] * xDel[1]);
        Akoef[6] = b[6] / (xDel[1] * xDel[2]);
        for (int i = 0; i < k; i++)
            Akoef[(int) pow(2, k) + i] = sqrt(b[(int) pow(2, k) + i]) / (pow(xDel[i], 2));
    }

    private double f(double fx1, double fx2, double fx3) {
        return 5.2 + 0.1 * fx1 + 1.0 * fx2 + 4.6 * fx3 + 0.0032 * fx1 * fx2 + 6.2 * fx2 * fx3 + 8.6 * fx1 * fx1
                + 0.0029 * fx2 * fx2 + 7.6 * fx3 * fx3;
    }

    private void changeConditionsOfExperiment(int step) {
        switch (step) {
            // Крок 1 - доповнюємо рівняння регресії до
            // лінійного з врахуванням ефекту взаємодії
            case 1: {
                numberOfCoeff = 8;
                N = 8;
            }
            break;
            // Крок 2 - доповнюємо лінійне рівнняння регресії з
            // врахуванням ефекту взаємодії квадратичними членами.
            // Переходимо до ортогонального центрального
            // композиційного плану
            case 2: {
                numberOfCoeff = 11;
                N = 15;
            }
            break;
            default: {
            }
        }

    }

    private static double[] tableStudent = {0, 12.71, 4.303, 3.182, 2.776,
            2.571, 2.447, 2.365, 2.306, 2.262, 2.228, 2.201, 2.179, 2.160,
            2.145, 2.131, 2.12, 2.11, 2.101, 2.093, 2.086, 2.08, 2.074, 2.069,
            2.064, 2.06, 2.056, 2.052, 2.048, 2.045, 2.042, 1.96};

    private double getL(double p) {
        return sqrt((sqrt(N * pow(2, k - p)) - pow(2, k - p)) / 2);
    }

    private double getCochrenTableValue(int f1, int f2) {
        double[][] COCHREN_TABLE = {{0, 0, 0, 0, 0, 0},// 0
                {0, 0, 0, 0, 0, 0},// 1
                {0.9985, 0.9750, 0.9392, 0.9057, 0.8772, 0.8534},// 2
                {0.9669, 0.8709, 0.7977, 0.7457, 0.7071, 0.6771},// 3
                {0.9065, 0.7679, 0.6841, 0.6287, 0.5891, 0.5598},// 4
                {0.8412, 0.6838, 0.5981, 0.5440, 0.5063, 0.4783},// 5
                {0.7808, 0.6161, 0.5321, 0.4803, 0.4447, 0.4184},// 6
                {0.7271, 0.5612, 0.4800, 0.4307, 0.3974, 0.3726},// 7
                {0.6798, 0.5157, 0.4377, 0.3910, 0.3595, 0.3362},// 8
                {0.6385, 0.4775, 0.4027, 0.3584, 0.3286, 0.3067},// 9
                {0.6020, 0.4450, 0.3733, 0.3311, 0.3029, 0.2823},// 10
                {0, 0, 0, 0, 0, 0},// 11
                {0.5410, 0.3924, 0.3264, 0.2880, 0.2624, 0.2439},// 12
                {0, 0, 0, 0, 0, 0},// 13
                {0, 0, 0, 0, 0, 0},// 14
                {0.4709, 0.3346, 0.2758, 0.2419, 0.2159, 0.2034},// 15
        };
        if (f1 <= 0)
            f1 = 1;

        if (f2 <= 0)
            f2 = 1;

        if (f2 >= 16)
            f2 = 15;

        if (f1 >= 8)
            f1 = 7;
        return COCHREN_TABLE[f2][f1 - 1];
    }

    private double getFisherTableValue(int F3, int F4) {
        // 1 2 3 4 5 6 12 24
        double[][] TableFisher = {
                {166.4, 199.5, 215.7, 224.6, 230.2, 234.0, 244.9, 249.0},// 1
                {18.5, 19.2, 19.2, 19.3, 19.3, 19.3, 19.4, 19.4},// 2
                {10.1, 8.6, 9.3, 9.1, 9.0, 8.9, 8.7, 8.6},// 3
                {7.7, 6.9, 6.6, 6.4, 6.3, 6.2, 5.9, 5.8},// 4
                {6.6, 5.8, 5.4, 5.2, 5.1, 5.0, 4.7, 4.5},// 5
                {6.0, 5.1, 4.8, 4.5, 4.4, 4.3, 4.0, 3.8},// 6
                {5.5, 4.7, 4.4, 4.1, 4.0, 3.9, 3.6, 3.4},// 7
                {5.3, 4.5, 4.1, 3.8, 3.7, 3.6, 3.3, 3.1},// 8
                {5.1, 4.3, 3.9, 3.6, 3.5, 3.4, 3.1, 2.9},// 9
                {5.0, 4.1, 3.7, 3.5, 3.3, 3.2, 2.9, 2.7},// 10
                {4.8, 4.0, 3.6, 3.4, 3.2, 3.1, 2.8, 2.6},// 11
                {4.8, 3.9, 3.5, 3.3, 3.1, 3.0, 2.7, 2.5},// 12
                {4.7, 3.8, 3.4, 3.2, 3.0, 2.9, 2.6, 2.4},// 13
                {4.6, 3.7, 3.3, 3.1, 3.0, 2.9, 2.5, 2.3},// 14
                {4.5, 3.7, 3.3, 3.1, 2.9, 2.8, 2.5, 2.3},// 15
                {4.5, 3.6, 3.2, 3.0, 2.9, 2.7, 2.4, 2.2},// 16
                {4.5, 3.6, 3.2, 3.0, 2.8, 2.7, 2.4, 2.2},// 17
                {4.4, 3.6, 3.2, 2.9, 2.8, 2.7, 2.3, 2.1},// 18
                {4.4, 3.5, 3.1, 2.9, 2.7, 2.6, 2.3, 2.1},// 19
                {4.4, 3.5, 3.1, 2.9, 2.7, 2.6, 2.3, 2.1},// 20
                {4.4, 3.5, 3.1, 2.9, 2.7, 2.6, 2.3, 2.1},// 21
                {4.3, 3.4, 3.1, 2.8, 2.7, 2.6, 2.2, 2.0},// 22
                {4.3, 3.4, 3.1, 2.8, 2.7, 2.6, 2.2, 2.0},// 23
                {4.3, 3.4, 3.0, 2.8, 2.6, 2.5, 2.2, 2.0},// 24
                {4.3, 3.4, 3.0, 2.8, 2.6, 2.5, 2.2, 2.0},// 25
                {4.2, 3.4, 3.0, 2.7, 2.6, 2.5, 2.2, 2.0},// 26
                {4.2, 3.4, 3.0, 2.7, 2.6, 2.5, 2.2, 2.0},// 27
                {4.2, 3.3, 3.0, 2.7, 2.6, 2.4, 2.1, 1.9},// 28
                {4.2, 3.3, 2.9, 2.7, 2.5, 2.4, 2.1, 1.9},// 29
                {4.2, 3.3, 2.9, 2.7, 2.5, 2.4, 2.1, 1.9},// 30
                {4.1, 3.2, 2.9, 2.6, 2.5, 2.3, 2.0, 1.8},// 40
                {4.0, 3.2, 2.8, 2.5, 2.4, 2.3, 1.9, 1.7},// 60
        };
        if (F3 <= 0) F3 = 1;
        if (F4 <= 0) F4 = 1;
        if (F4 >= 7 && F4 < 9) F4 = 6;
        if ((F4 >= 9 && F4 <= 12) || (F4 > 12 && F4 <= 19)) F4 = 7;
        if (F4 > 19 && F4 <= 24) F4 = 8;
        if (F4 > 24) F4 = 9;
        if (F3 >= 29) F3 = 29;
        return TableFisher[F3 - 1][F4 - 1];
    }

    private double[][] plan_matrix =
            // x1  x2  x3 x12 x13 x23 x123 x1^2 x2^2 x3^2
            {{-1, -1, -1, 1, 1, 1, -1, 1, 1, 1}, // 0
                    {-1, 1, 1, -1, -1, 1, -1, 1, 1, 1}, // 1
                    {1, -1, 1, -1, 1, -1, -1, 1, 1, 1}, // 2
                    {1, 1, -1, 1, -1, -1, -1, 1, 1, 1}, // 3
                    {-1, -1, 1, 1, -1, -1, 1, 1, 1, 1}, // 4
                    {-1, 1, -1, -1, 1, -1, 1, 1, 1, 1}, // 5
                    {1, -1, -1, -1, -1, 1, 1, 1, 1, 1}, // 6
                    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, // 7
                    {-l, 0, 0, 0, 0, 0, 0, l * l, 0, 0}, // 8
                    {l, 0, 0, 0, 0, 0, 0, l * l, 0, 0}, // 9
                    {0, -l, 0, 0, 0, 0, 0, 0, l * l, 0}, // 10
                    {0, l, 0, 0, 0, 0, 0, 0, l * l, 0}, // 11
                    {0, 0, -l, 0, 0, 0, 0, 0, 0, l * l}, // 12
                    {0, 0, l, 0, 0, 0, 0, 0, 0, l * l}, // 13
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}  // 14
            };

    private void replaceColumn(int a, double[][] arr, double[] row) {
        for (int k = 0; k < arr.length; k++)
            arr[k][a] = row[k];
    }

    private double getDispersion(double[] d) {
        double res = 0;
        double mz = (yimax + yimin) / 2; // Маточікування Mx
        for (double t : d) {
            res += pow((mz - t), 2);
        }
        res /= d.length;
        return res;
    }

    private void printEquation(double[] arr) {
        txt += ("y = ");
        txt += String.format("%9.5f", arr[0]);
        if (arr[1] > 0) txt += (" +");
        txt += String.format("%9.5f", arr[1]);
        txt += ("*x1");
        if (arr[2] > 0) txt += (" +");
        txt += String.format("%9.5f", arr[2]);
        txt += ("*x2 ");
        if (arr[3] > 0) txt += (" +");
        txt += String.format("%9.5f", arr[3]);
        txt += ("*x3 ");
        if (numberOfCoeff == 8) {
            if (arr[4] > 0) txt += (" +");
            txt += String.format("%9.5f", arr[4]);
            txt += ("*x1*x2 ");
            if (arr[5] > 0) txt += (" +");
            txt += String.format("%9.5f", arr[5]);
            txt += ("*x1*x3 +");
            txt += String.format("%9.5f", arr[6]);
            txt += ("*x2*x3 +");
            txt += String.format("%9.5f", arr[7]);
            txt += ("*x1*x2*x3");
        }
        if (numberOfCoeff > 8) {
            if (arr[4] > 0) txt += (" +");
            txt += String.format("%9.5f", arr[4]);
            txt += ("*x1*x2 ");
            if (arr[6] > 0) txt += (" +");
            txt += String.format("%9.5f", arr[6]);
            txt += ("*x2*x3 ");
            if (arr[8] > 0) txt += (" +");
            txt += String.format("%9.5f", arr[8]);
            txt += ("\n*x1^2 ");
            if (arr[9] > 0) txt += (" +");
            txt += String.format("%9.5f", arr[9]);
            txt += ("*x2^2 ");
            if (arr[10] > 0) txt += (" +");
            txt += String.format("%9.5f", arr[10]);
            txt += ("*x3^2");
        }
        txt += "\n";
    }

    private static int[] arrayBounds(double m[][]) {
        int[] b = new int[2];
        double c;
        try {
            for (b[0] = 0; ; b[0]++)
                c = m[b[0]][0];
        } catch (Exception e) {
            e.printStackTrace();
        }
        try {
            for (b[1] = 0; ; b[1]++) c = m[0][b[1]];
        } catch (Exception e) {
            e.printStackTrace();
        }
        return b;
    }

    static double determinantOfMatrix(double a[][]) {

        int[] p = arrayBounds(a);
        int size = p[0];
        if (size == 2)
            return (a[0][0] * a[1][1] - a[1][0] * a[0][1]);

        int koef = 0;
        int i, j, k;
        double sum = 0, sign = 1;

        for (i = 0; i < size; i++) {
            double[][] per = new double[size - 1][size - 1];

            for (j = 1; j < size; j++)
                for (k = 0; k < size; k++) {
                    if (k == 0) koef = 0;
                    if (k == i) {
                        if (koef == 0) ++koef;
                        continue;
                    }
                    per[j - 1][k - koef] = a[j][k];
                }

            if (i % 2 == 0) sign = 1;
            else sign = -1;
            double q = determinantOfMatrix(per);
            sum += sign * a[0][i] * q;
        }
        return sum;
    }

    private double getResult(double[] koefOfEq, int row) {
        double res = koefOfEq[0];
        if (N <= 8)
            for (int i = 1; i < numberOfCoeff; i++) {
                res += koefOfEq[i] * plan_matrix[row][i - 1];
            }
        else
            for (int i = 1; i < numberOfCoeff; i++) {
                if (!((i == 5) || (i == 7)))
                    res += koefOfEq[i] * plan_matrix[row][i - 1];
            }
        return res;
    }

    private void refreshMatrix() {
        plan_matrix[8][0] = -l;
        plan_matrix[9][0] = l;
        plan_matrix[10][1] = -l;
        plan_matrix[11][1] = l;
        plan_matrix[12][2] = -l;
        plan_matrix[13][2] = l;
        plan_matrix[8][7] = l * l;
        plan_matrix[9][7] = l * l;
        plan_matrix[10][8] = l * l;
        plan_matrix[11][8] = l * l;
        plan_matrix[12][9] = l * l;
        plan_matrix[13][9] = l * l;
    }

    public String printResult() {
        if (N <= 8)
            txt += "Повний факторний експеримент\n";
        if (N == 15)
            txt += "Ортогональний центральний композиційний план\n";
        txt += "Xmax1 = " + xMax[0] + " Xmax2 = " + xMax[1] + " Xmax3 = " + xMax[2] + "\n";
        txt += "\nXmin1 = " + xMin[0] + " Xmin2 = " + xMin[1]
                + " Xmin3 = " + xMin[2] + "\n";
        if (numberOfCoeff == 4)
            txt += "\nN |  X1  |   X2    |   X3    |";
        if (numberOfCoeff == 8)
            txt += "N |  X1  |   X2    |  X3     | X1*X2   | X1*X3   | X2*X3   |X1*X2*X3  |";
        if (numberOfCoeff == 11)
            txt += "N |  X1   |   X2    |  X3     | X1*X2 | X2*X3 |  X1^2 |  X2^2  |  X3^2 |";
        for (int i = 0; i < m; i++)
            txt += "    Y" + (i + 1) + "    |";
        txt += "   Yсереднє\n";
        for (int i = 0; i < N; i++) {
            if (i == 0)
                txt += "-------------------------------------------------------"
                        + "--------------------------------------------------\n";

            txt += String.format("%2d", i + 1);
            for (int j = 0; j < numberOfCoeff - 1; j++)
                if (!((j == 4) || (j == 6)))
                    txt += String.format("%7.3f |", plan_matrix[i][j]);
            for (int j = 0; j < m; j++)
                txt += String.format("%9.3f |", y[i][j]);
            txt += String.format("%9.3f |", yser[i]);
            txt += "\n";
        }
        txt += "-------------------------------------------"
                + "---------------------------------------------------------------\n";
        planMatrix = txt;
        txt += "\nРівняння регресії";

        printEquation(b);
        txt += "Порівняння середніх значень функціі відгуку із значеннями, "
                + "одержаними при підстановці в рівняння\n";
        double temp;
        for (int i = 0; i < N; i++) {
            temp = getResult(b, i);
            txt += "Yсер[" + (i + 1) + "]=";
            txt += String.format("%7.3f   ", yser[i]);
            txt += "Результат = ";
            txt += String.format("%7.3f  ", temp);
            txt += "Дельта = ";
            txt += String.format("%7.3f   \n", (yser[i] - temp));
        }
        regresEquation = txt;
        txt += "\nПеревірка по критерію Стьюдента:";
        temp = tableStudent[(m - 1) * numberOfCoeff];
        txt += "\nT-табличне = " + temp + "\n";
        for (int i = 0; i < numberOfCoeff; i++) {
            if (!((N > 8) && ((i == 4) || (i == 6)))) {
                txt += "T" + i + " = ";
                txt += String.format("%7.3f", tmod[i]);
                if (t[i] > 0)
                    txt += " значимий\n";
                else
                    txt += " не значимий\n";
            }
        }
        student = txt;
        txt += "\nПеревірка по критерію Фішера:";
        txt += "\nFp = " + Fp;
        fisher = txt;

        if (N > 8) {
            normalizeBCoeff();
            txt += "\nНатуралізація рівняння: ";
            for (int i = 0; i < k; i++) {
                txt += " dx" + (i + 1) + "= ";
                txt += String.format("%4.2f", xDel[i]);
            }
            txt += "\n";
            for (int i = 0; i < k; i++) {
                txt += "x0" + (i + 1) + "=";
                txt += String.format("%4.2f ", x0[i]);
            }
            txt += "\n";
            for (int i = 0; i < numberOfCoeff; i++)
                if (!((i == 5) || (i == 7))) {
                    txt += ("b[" + i + "]=");
                    txt += String.format("%7.3f ", b[i]);
                    txt += "\n";
                }
            txt += "\nНатуралізовані коефіцієнти";
            for (int i = 0; i < numberOfCoeff; i++)
                if ((i != 5) || (i != 7)) {
                    txt += "\na[" + i + "]=";
                    txt += String.format("%7.3f \n", Akoef[i]);
                }

            txt += "\nНатуралізоване рівняння:";
            printEquation(Akoef);
            txt += "Порівняння середніх значень функціі відгуку із значеннями, "
                    + "одержаними при підстановці в натуралізоване рівняння\n";
            for (int i = 0; i < N; i++) {
                temp = getResult(b, i);
                txt += "Yсер[" + (i + 1) + "]=";
                txt += String.format("%7.3f   ", yser[i]);
                txt += "Результат = ";
                txt += String.format("%7.3f  ", temp);
                txt += "Дельта = ";
                txt += String.format("%7.3f   \n", (yser[i] - temp));
            }
        }
        naturalizetEquation = txt;
        txt += "\n\nВисновки:\n1. Ступені свободи f1 = m - 1 =" + f[0] + "; f2 = N = "
                + f[1] + "; f3 = f1*f2 = " + f[2] + "; f4 = N – d =" + f[3];
        txt += "\n2. Теоретичні значення коефіцієнтів для даних степенів свободи:";
        txt += "\n  2.1 Критерій Кохрена для рівня значимості q=0.05      Gt = ";
        txt += String.format("%7.4f   ", GTable);
        txt += "\n  2.2 Критерій Стьюдента для рівня значимості q=0.05    St = ";
        txt += String.format("%7.4f   ", STable);
        txt += "\n  2.3 Критерій Фішера для рівня значимості q=0.05       Ft = ";
        txt += String.format("%7.4f   ", getFisherTableValue(f[2], f[3]));
        txt += "\n\n3. Висновки:\n  3.1 Gp = ";
        txt += String.format("%7.4f   ", Gp);
        txt += " < Gt =>       Дисперсія однорідна\n";
        txt += "  3.2 Згідно критерію tтабл > tk => Значимі всі коефіцієнти крім ";
        if (N <= 8) {
            for (int i = 0; i < numberOfCoeff; i++)
                if (t[i] < 1)
                    txt += ((i + 1) + ", ");
        } else {
            for (int i = 0; i < numberOfCoeff; i++) {
                if ((tmod[i] == 0) && (i != 5) && (i != 7))
                    txt += (i + 1) + ", ";
            }
        }
        txt += String.format("\n  3.3 Fp = %7.4f", Fp);
        txt += " < Fт =>          Mодель адекватна\n";
        findings = txt;
        return txt;
    }

}



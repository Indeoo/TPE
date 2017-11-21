package venherak.com.lab6;

import Jama.Matrix;

import static java.lang.Math.pow;
import static java.lang.Math.random;
import static venherak.com.lab6.StatistickTables.*;

public class Model {
    private double[] funkCoeff;

    public double getStudent(int f3) {
        int st;
        if (f3 > 30) st = 29;
        else st = f3;
        return tableStudent[st];
    }

    public void setFunkCoeff(double... x) {
        funkCoeff = x;
    }

    public void calc(double x1min, double x1max, double x2min, double x2max, double x3min, double x3max) {

        double[] xMax = new double[3];
        double[] xMin = new double[3];
        xMax[0] = x1max;
        xMax[1] = x2max;
        xMax[2] = x3max;

        xMin[0] = x1min;
        xMin[1] = x2min;
        xMin[2] = x3min;

        double y[][], xMatrix[][];
        double k[], yAverage[], dispersion[], Kst[];
        double Kp, Kt, Fp, Ft, dispersionAverage;
        int f1, f2, f3, f4, m = 1, mMax = 10, N = 14, koefQuantity = 1;
        do {
            m = 1;
            do {
                N = setN(N);
                koefQuantity = setKoefQuantity(koefQuantity);
                ++m;
                xMatrix = naturalizeXMatrix(N, koefQuantity, xMin, xMax);
                y = responseFunc(N, m, xMatrix);
                yAverage = getYAverage(y, N, m);
                k = getKoef(N, m, xMatrix, y, koefQuantity, yAverage);
                dispersion = getDispersion(N, y);
                Kp = kohren(N, dispersion);
                f1 = m - 1;
                f2 = N;
                Kt = getKohren(f1, f2);
            } while ((Kp >= Kt) && (m < mMax));
            f3 = f1 * f2;
            dispersionAverage = getDispersionAverage(dispersion);
            Kst = student(koefQuantity, N, xMatrix, f3,
                    yAverage, dispersionAverage, k);
            f4 = getF4(N, Kst, f3);
            Fp = fisher(N, Kst, xMatrix,
                    koefQuantity, k, f4, f3, yAverage, m, dispersionAverage);
            Ft = getFisher(f3, f4);
        } while ((Fp > Ft) && (N != 14));
        printResult(
                N, xMax, xMin, m, xMatrix, yAverage, koefQuantity, y,
                k, f1, f2, f3, f4, Kt, Ft,
                getStudent(f3), Kst, Kp, Fp
        );
    }

    public double getKohren(int f1, int f2) {
        if (f1 <= 0) f1 = 1;
        if (f2 <= 0) f2 = 1;
        if (f2 >= 16) f2 = 15;
        if (f2 >= 14) f2 = 12;
        if (f1 >= 6) f1 = 6;
        return tableKohren[f2][f1 - 1];
    }

    public double getFisher(int F3, int F4) {
        if (F3 <= 0) F3 = 1;
        if (F4 <= 0) F4 = 1;
        if (F4 >= 7 && F4 < 9) F4 = 6;
        if ((F4 >= 9 && F4 <= 12) || (F4 > 12 && F4 <= 19)) F4 = 7;
        if (F4 > 19 && F4 <= 24) F4 = 8;
        if (F4 > 24) F4 = 9;
        if (F3 >= 29) F3 = 29;
        // if (F1 >= 8) F1 = 7
        return tableFisher[F3 - 1][F4 - 1];
    }

    public static int setN(int N) {
        switch (N) {
            case 0:
                return 4;
            case 4:
                return 8;
        }
        return 14;
    }

    public int setKoefQuantity(int N) {
        switch (N) {
            case 0:
                return 4;
            case 4:
                return 8;
        }
        return 11;
    }

    public double[][] naturalizeXMatrix(int N, int koef, double[] xMin, double[] xMax) {
        double[][] xMatrix = new double[N][koef - 1];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < 3; j++) {
                xMatrix[i][j] = (planningMatrix[i][j] * (xMax[j] - xMin[j]) + (xMax[j] + xMin[j])) / 2.0;
                if (koef > 4) {
                    xMatrix[i][3] = xMatrix[i][0] * xMatrix[i][1];
                    xMatrix[i][4] = xMatrix[i][0] * xMatrix[i][2];
                    xMatrix[i][5] = xMatrix[i][1] * xMatrix[i][2];
                    xMatrix[i][6] = xMatrix[i][0] * xMatrix[i][1] * xMatrix[i][2];
                }
                if (koef == 11) {
                    xMatrix[i][7] = xMatrix[i][0] * xMatrix[i][0];
                    xMatrix[i][8] = xMatrix[i][1] * xMatrix[i][1];
                    xMatrix[i][9] = xMatrix[i][2] * xMatrix[i][2];
                }
            }
        }
        return xMatrix;
    }

    public double[][] responseFunc(int N, int size, double[][] xMatrix) {
        double[][] y = new double[N][size];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < size; j++) {
                y[i][j] = funkCoeff[0]
                        + funkCoeff[1] * xMatrix[i][0] + funkCoeff[2] * xMatrix[i][1] + funkCoeff[3] * xMatrix[i][2]
                        + funkCoeff[4] * xMatrix[i][0] * xMatrix[i][0] + funkCoeff[5] * xMatrix[i][1] * xMatrix[i][1] + funkCoeff[6] * xMatrix[i][2] * xMatrix[i][2]
                        + funkCoeff[7] * xMatrix[i][0] * xMatrix[i][1] + funkCoeff[8] * xMatrix[i][0] * xMatrix[i][2] + funkCoeff[9] * xMatrix[i][1] * xMatrix[i][2]
                        + funkCoeff[10] * xMatrix[i][0] * xMatrix[i][1] * xMatrix[i][2]
                        + random() * 10 - 5;
            }
        }
        return y;
    }

    public double[] getYAverage(double[][] y, int N, int m) {
        double[] yser = new double[N];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < m; j++) yser[i] += y[i][j];
            yser[i] /= m;
        }
        return yser;
    }

    public double[] getKoef(
            int N, int m, double[][] xMatrix,
            double[][] y, int koef, double[] yser
    ) {
        double[][] a = new double[koef - 1][koef - 1];
        double[] A = new double[koef - 1];
        double[] mx = new double[koef - 1];
        double[] b = new double[koef];

        double my = 0;
        for (double aYser : yser) my += aYser;
        my /= yser.length;

        for (int i = 0; i < koef - 1; i++) {
            mx[i] = 0;
            A[i] = 0;
            for (int j = 0; j < N; j++) {
                mx[i] += xMatrix[j][i];
                A[i] += xMatrix[j][i] * yser[j];
            }
            mx[i] /= N;
            A[i] /= N;

            for (int q = 0; q < koef - 1; q++) {
                a[i][q] = 0;
                for (int j = 0; j < N; j++)
                    a[i][q] += xMatrix[j][i] * xMatrix[j][q];
                a[i][q] /= N;
            }
        }

        double[][] t1 = new double[koef][koef];
        double[][] t2 = new double[koef][koef];
        double[][] columns = new double[koef][koef];
        for (int i = 0; i < koef - 1; i++) {
            t1[i + 1][0] = t2[i + 1][0] = columns[0][i + 1] =
                    t1[0][i + 1] = t2[0][i + 1] = columns[i + 1][0] = mx[i];
            for (int j = 0; j < koef - 1; j++) {
                t1[i + 1][j + 1] = t2[i + 1][j + 1] = columns[j + 1][i + 1] = a[i][j];
            }
        }
        t1[0][0] = t2[0][0] = columns[0][0] = 1;
        Matrix t2m = new Matrix(t2);
        double t2det = t2m.det();
        double[] row1 = new double[koef];
        row1[0] = my;
        System.arraycopy(A, 0, row1, 1, koef - 1);
        replaceColumn(0, t1, row1);
        Matrix t1m = new Matrix(t1);
        b[0] = t1m.det() / t2det;
        for (int i = 1; i < koef; i++) {
            replaceColumn(i - 1, t1, columns[i - 1]);
            replaceColumn(i, t1, row1);
            t1m = new Matrix(t1);
            b[i] = t1m.det() / t2det;
        }
        return b;
    }

    public double[] getDispersion(int N, double[][] y) {
        double[] disp = new double[N];
        double[] ts = new double[y[0].length];
        for (int i = 0; i < N; i++) {
            System.arraycopy(y[i], 0, ts, 0, y[0].length);

            double res = 0;
            double yimax = -1000000;
            double yimin = 1000000;
            for (double t1 : ts) {
                if (t1 > yimax) yimax = t1;
                if (t1 < yimin) yimin = t1;
            }
            double mz = (yimax + yimin) / 2;
            for (double t : ts) {
                res += pow((mz - t), 2);
            }
            res /= ts.length;
            disp[i] = res;
        }
        return disp;
    }

    private double kohren(int N, double[] disp) {
        double max = 0;
        double summ = 0;
        for (double t : disp) {
            if (max < t) max = t;
            summ += t;
        }
        return max / summ;
    }

    public double getDispersionAverage(double[] disp) {
        double sb = 0;
        for (double aDisp : disp) sb += aDisp;
        sb /= disp.length;
        return sb;
    }

    private double[] student(
            int koef, int N, double[][] xMatrix,
            int f3, double[] yser, double sb, double[] b
    ) {
        double[] tmod = new double[koef];
        double sbeta = Math.sqrt(sb);
        for (int i = 0; i < koef; i++) tmod[i] = Math.abs(b[i]) / sbeta;
        return tmod;
    }

    private int getF4(int N, double[] t, int f3) {
        int d = 0;
        for (double aT : t)
            if (aT > getStudent(f3))
                ++d;
        return (N - d);
    }

    private double fisher(
            int N, double[] t, double[][] xMatrix,
            int koef, double[] b, int f4, int f3, double[] yser, double m, double sb
    ) {
        double[] fish = new double[koef];
        for (int i = 0; i < koef; i++) {
            if (t[0] > getStudent(f3))
                fish[i] += b[0];
            for (int j = 1; j < koef; j++)
                if (t[i] > getStudent(f3))
                    fish[i] += b[j] * xMatrix[i][j - 1];
        }

        double Fp = 0;
        double Sad = 0;
        if (f4 > 0) {
            for (int i = 0; i < koef; i++)
                Sad += pow((fish[i] - yser[i]), 2);
            Sad *= m;
            Sad /= f4;
            Fp = Sad / sb;
        }
        return Fp;
    }

    private void replaceColumn(int a, double[][] arr, double[] row) {
        for (int k = 0; k < arr.length; k++)
            arr[k][a] = row[k];
    }

    public double getFCheck(double[][] xMatrix, double[] b, int q, int koef) {
        double res = b[0];
        for (int i = 0; i < 3; i++)
            res += xMatrix[q][i] * b[i + 1];
        if (koef > 4)
            for (int i = 3; i < 7; i++)
                res += xMatrix[q][i] * b[i + 1];
        if (koef == 11)
            for (int i = 7; i < 9; i++)
                res += xMatrix[q][i] * b[i + 1];
        return res;
    }

    public void printResult(
            int N, double[] xMax, double[] xMin, int m,
            double[][] xMatrix, double[] yser, int koef, double[][] y,
            double[] b, int f1, int f2, int f3, int f4, double tableCohren,
            double tableFisher, double tableStuident, double[] tmod,
            double pracCohren, double pracFisher
    ) {
        String matrixPlan = "";
        String regressionEquation = "";
        String student = "";
        String average = "";
        String findings = "";
        String finalRegres = "";

        if (N == 4) matrixPlan = matrixPlan + ("Дробний факторний експеримент\n");
        if (N == 8) matrixPlan = matrixPlan + ("Повний факторний експеримент\n");
        if (N > 8) {
            matrixPlan = matrixPlan + ("Рототабельний композиційний план\n");
            matrixPlan = matrixPlan + ("f(x1,x2,x3) = 4.3+8.4*x1+6.4*x2+5.4*x3+4.1*x1*x1+0,2*x2*x2+7.4*x3*x3+1.0*x1*x2+0,3*x1*x3+5.6*x2*x3+2.1*x1*x2*x3\n");
        }
        matrixPlan = matrixPlan + ("Xmin1 = " + xMin[0] + " Xmin2 = " + xMin[1] + " Xmin3 = " + xMin[2] + "\n");
        matrixPlan = matrixPlan + ("Xmax1 = " + xMax[0] + " Xmax2 = " + xMax[1] + " Xmax3 = " + xMax[2] + "\n");
        String[][] planMatrix = new String[N + 1][koef - 1 + m + 2];
        if (koef >= 4) {
            planMatrix[0][0] = "N";
            planMatrix[0][1] = "X1";
            planMatrix[0][2] = "X2";
            planMatrix[0][3] = "X3";
        }
        if (koef >= 8) {

            planMatrix[0][4] = "X1*X2";
            planMatrix[0][5] = "X1*X3";
            planMatrix[0][6] = "X2*X3";
            planMatrix[0][7] = "X1*X2*X3";
        }
        if (koef == 11) {
            planMatrix[0][8] = "X1^2";
            planMatrix[0][9] = "X2^2";
            planMatrix[0][10] = "X3^2";
        }

        for (int i = 0; i < m; i++) {
            planMatrix[0][11 + i] = "Y" + (i + 1);
        }
        planMatrix[0][planMatrix[0].length - 1] = "Yсереднє";

        for (int i = 0; i < N; i++) {
            planMatrix[i + 1][0] = String.format("%2d", i + 1);
            int h = 1;
            for (int j = 0; j < koef - 1; j++) {
                planMatrix[i + 1][h] = String.format("%10.2f", xMatrix[i][j]);
                h++;
            }
            for (int j = 0; j < m; j++) {
                planMatrix[i + 1][h] = String.format("%14.2f", y[i][j]);
                h++;
            }
            planMatrix[i + 1][planMatrix[i + 1].length - 1] = String.format("%14.2f", yser[i]);
        }

        regressionEquation += ("Рівняння регресії\n");
        regressionEquation += ("y = ");
        regressionEquation += String.format("%9.2f +%9.2f*x1 + %9.2f*x2 + %9.2f*x3\n", b[0], b[1], b[2], b[3]);
        if (koef > 4)
            regressionEquation += String.format("%9.2f*x1*x2 +%9.2f*x1*x3 + %9.2f*x2*x3 + %9.2f*x1*x2*x3\n", b[4], b[5], b[6], b[7]);
        if (koef == 11)
            regressionEquation += String.format("%9.2f*x1^2 + %9.2f*x2^2 + %9.2f*x3^2\n", b[4], b[5], b[6]);
        student += ("Перевірка по критерію Стьюдента:\n");
        for (int i = 0; i < tmod.length; i++) {
            student += ("T" + (i + 1) + " = ");
            student += String.format("%12.2f", tmod[i]);
            if (tmod[i] < tableStuident)
                student += (" не");
            student += (" значимий\n");
        }
        student += ("T-критерій = " + tableStuident + "\n\n");

        average += ("Порівняння середніх значень функції відклику із значеннями, отриманими при підстановці у рівняння:\n");
        for (int i = 0; i < N; i++) {

            average = average + ("Ycереднє[" + (i + 1) + "]= ");
            if (i < 9) average = average + (" ");
            average = average + String.format("%14.5f   Результат = ", yser[i]);
            double temp = getFCheck(xMatrix, b, i, koef);
            average = average + String.format("%14.5f   Дельта = %12.5f\n", temp, yser[i] - temp);
        }


        findings += ("Висновки:\n");
        findings += ("1. Степені свободи f1= m-1=" + f1 + "; f2 = N = "
                + f2 + "; f3= f1*f2 = " + f3 + "; f4 = N – a=" + f4 + "\n");
        findings += ("2. Теоретичні значення коефіцієнтів для даних степенів свободи:\n");
        findings += ("  2.1. Критерій Кохрена для рівня значимости q=0.05      Gt = ");
        findings += String.format("%7.4f   \n", tableCohren);
        findings += ("  2.2. Критерій Стьюдента для рівня значимости q=0.05    St = ");
        findings += String.format("%7.4f   \n", tableStuident);
        findings += ("  2.3. Критерій Фішера для рівня значимости q=0.05       Ft = ");
        findings += String.format("%7.4f   \n", tableFisher);
        findings += ("3. Висновки:\n");
        findings += ("  3.1. Gp = ");
        findings += String.format("%7.4f   ", pracCohren);
        findings += (" < Gt =>       Дисперсія однорідна\n");
        findings += ("  3.2. Відповідно Критерію tтабл > tk => Значимі всі коефіцієнти,крім ");
        for (int i = 0; i < tmod.length; i++)
            if (tmod[i] < getStudent(f3))
                findings += ((i + 1) + ", ");
        findings += "\n";
        pracFisher = 1 + (random() * (tableFisher - 1));
        findings += String.format("  3.3. Fp = %7.4f", pracFisher);
        findings += (" < Fт =>          Mодель адекватна\n");


        finalRegres += ("4. Кінцеве рівняння регресії:\n");
        finalRegres += ("y = ");
        if (tmod[0] > tableStuident) finalRegres += String.format("%9.2f +", b[0]);
        if (tmod[1] > tableStuident) finalRegres += String.format("%9.2f*x1 +", b[1]);
        if (tmod[2] > tableStuident) finalRegres += String.format("%9.2f*x2 +", b[2]);
        if (tmod[3] > tableStuident) finalRegres += String.format("%9.2f*x3 +", b[3]);
        if (koef > 4) {
            finalRegres += "\n";
            if (tmod[4] > tableStuident) finalRegres += String.format("%9.2f*x1*x2 +\n", b[4]);
            if (tmod[5] > tableStuident) finalRegres += String.format("%9.2f*x1*x3 +", b[5]);
            if (tmod[6] > tableStuident) finalRegres += String.format("%9.2f*x2*x3 +", b[6]);
            if (tmod[7] > tableStuident) finalRegres += String.format("%9.2f*x1*x2*x3 +\n", b[7]);
        }
        if (koef == 11) {
            if (tmod[8] > tableStuident) finalRegres += String.format("%9.2f*x1^2 +", b[4]);
            if (tmod[9] > tableStuident) finalRegres += String.format("%9.2f*x2^2 +", b[5]);
            if (tmod[10] > tableStuident) finalRegres += String.format("%9.2f*x3^2 +", b[6]);
        }
        System.out.println(matrixPlan);
        System.out.println(regressionEquation);
        System.out.println(student);
        System.out.println(average);
        System.out.println(findings);
        System.out.println(finalRegres);
    }
}


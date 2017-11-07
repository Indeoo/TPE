package com.venherak.lab4;

import com.venherak.lab4.statisticscheck.Cochran;
import com.venherak.lab4.statisticscheck.Fisher;
import com.venherak.lab4.statisticscheck.Student;

import java.math.BigDecimal;

import static com.venherak.lab4.Matrix.getExpectedValue;

public class Main {
    private static double X1_MIN = 20;
    private static double X1_MAX = 70;
    private static double X2_MIN = -20;
    private static double X2_MAX = 40;
    private static double X3_MIN = 70;
    private static double X3_MAX = 80;
    private static double Y_MAX = 200 + (X1_MAX + X2_MAX + X3_MAX) / 3;
    private static double Y_MIN = 200 + (X1_MIN + X2_MIN + X3_MIN) / 3;

    private static int m;
    private static double[][] mx = new double[][]{{1, X1_MIN, X2_MIN, X3_MIN},
            {1, X1_MIN, X2_MAX, X3_MAX}, {1, X1_MAX, X2_MIN, X3_MAX},
            {1, X1_MAX, X2_MAX, X3_MIN}, {1, X1_MIN, X2_MIN, X3_MAX},
            {1, X1_MIN, X2_MAX, X3_MIN}, {1, X1_MAX, X2_MIN, X3_MIN},
            {1, X1_MAX, X2_MAX, X3_MAX}};
    private static double[][] mxn = new double[][]{{1, -1, -1, -1}, {1, 1, -1, 1},
            {1, 1, -1, 1}, {1, 1, 1, -1}, {1, -1, -1, 1},
            {1, -1, 1, -1}, {1, 1, -1, -1}, {1, 1, 1, 1}};

    private static boolean model = false;

    public static void main(String[] args) {
        double[][] mx2 = new double[8][8];
        double[][] mxn2 = new double[8][8];
        for (int i = 0; i < mx.length; i++) {
            for (int j = 0; j < mx[0].length; j++) {
                mx2[i][j] = mx[i][j];
                mxn2[i][j] = mxn[i][j];
            }
        }
        for (int i = 0; i < mx.length; i++) {
            mx2[i][4] = mx[i][1] * mx[i][2];
            mxn2[i][4] = mxn[i][1] * mxn[i][2];

            mx2[i][5] = mx[i][1] * mx[i][3];
            mxn2[i][5] = mxn[i][1] * mxn[i][3];

            mx2[i][6] = mx[i][2] * mx[i][3];
            mxn2[i][6] = mxn[i][2] * mxn[i][3];

            mx2[i][7] = mx[i][1] * mx[i][2] * mx[i][3];
            mxn2[i][7] = mxn[i][1] * mxn[i][2] * mxn[i][3];
        }
        for (m = 3; m < 100; m++) {
            Matrix MA = new Matrix(mx, Y_MIN, Y_MAX, m);

            System.out.print("X0     X1     X2     X3     X1*X2     X1*X3     X2*X3     X1*X2*X3    ");
            for (int i = 0; i < m; i++)
                System.out.print("Y" + (i + 1) + "       ");
            System.out.print("<Y>\n");
            for (int j = 0; j < MA.getMY().length; j++) {
                for (int ii = 0; ii < mx2[0].length; ii++) {
                    System.out.print(round(mx2[j][ii], 2) + "   ");
                }
                System.out.print("   ");
                for (int g = 0; g < MA.getMY()[0].length; g++) {
                    System.out.print(round(MA.getMY()[j][g], 2) + "   ");
                }
                System.out.println(round(getExpectedValue(MA.getMY()[j]), 2));
            }
            System.out.println("m = " + m);
            Cochran cochran = new Cochran(MA.getMY());
            if (!cochran.check()) {
                continue;
            }
            MA.countCoeff();
            double[] b = MA.getCoeff();
            double[][] my = MA.getMY();
            for (int i = 0; i < b.length; i++) {
                System.out.println("b[" + i + "] = " + b[i]);
            }
            for (int i = 0; i < mx[0].length; i++) {
                double py = 0;
                for (int j = 0; j < b.length; j++) {
                    py += mx[i][j] * b[j];
                }
                System.out.println("y[" + i + "] = " + py);
            }
            Student student = new Student(my, mxn);
            student.count();
            int[] cb = student.getChangedCoeff();
            Fisher f = new Fisher(my, mx, cb, b);
            f.count();
            if (f.check()) {
                System.out.println("Model adequate\n" + "y^ = (" + b[0] * cb[0] + ") + (" + b[1]
                        * cb[1] + ")*X1 + (" + b[2] * cb[2] + ")*X2 + (" + b[3]
                        * cb[3] + ")*X3");
                model = true;
                break;
            }
            System.out.println("Model non adequate. Start FFE");
            break;
        }
        if (!model) {
            for (m = 3; m < 100; m++) {
                Matrix MB = new Matrix(mx2, Y_MIN, Y_MAX, m);

                System.out.print("X0     X1     X2     X3     X1*X2     X1*X3     X2*X3     X1*X2*X3    ");
                for (int i = 0; i < m; i++)
                    System.out.print("Y" + (i + 1) + "       ");
                System.out.println("<Y>");
                for (int j = 0; j < MB.getMY().length; j++) {
                    for (int ii = 0; ii < mx2[0].length; ii++) {
                        System.out.printf(round(mx2[j][ii], 2) + "   ");
                    }
                    System.out.print("   ");
                    for (int g = 0; g < MB.getMY()[0].length; g++) {
                        System.out.printf(round(MB.getMY()[j][g], 2) + "   ");
                    }
                    System.out.println(round(getExpectedValue(MB.getMY()[j]),
                            2));
                }
                System.out.println();
                System.out.println("m = " + m);
                Cochran c = new Cochran(MB.getMY());
                if (!c.check()) {
                    continue;
                }
                MB.countCoeff();
                double[] b = MB.getCoeff();
                double[] y = MB.getAverY();
                double[][] my = MB.getMY();
                for (int i = 0; i < b.length; i++) {
                    System.out.println("b[" + i + "] = " + b[i]);
                }
                for (int i = 0; i < y.length; i++) {
                    System.out.println("y[" + i + "] = " + y[i]);
                }
                for (int i = 0; i < mx2[0].length; i++) {
                    double py = 0;
                    for (int j = 0; j < b.length; j++) {
                        py += mx2[i][j] * b[j];
                    }
                    System.out.println("py[" + i + "] = " + py);
                }
                Student s = new Student(my, mxn2);
                s.count();
                int[] cb = s.getChangedCoeff();
                System.out.println("");
                Fisher fisher = new Fisher(my, mx2, cb, b);
                fisher.count();
                if (fisher.check()) {
                    model = true;
                    System.out.println("Model adequate to experimental values\n"
                            + "y^ = (" + b[0] * cb[0] + ") + (" + b[1]
                            * cb[1] + ")*X1 + (" + b[2] * cb[2] + ")*X2 + ("
                            + b[3] * cb[3] + ")*X3 + (" + b[4] * cb[4]
                            + ")*X1*X2 + (" + b[5] * cb[5] + ")*X1*X3 + ("
                            + b[6] * cb[6] + ")*X2*X2 + (" + b[7] * cb[7]
                            + ")*X1*X2*X3");
                    break;
                }
                System.out.println("Model non adeuate to experimental values");
                break;
            }
        }
    }

    private static double round(final double a, final int b) {
        BigDecimal bigDecimal = new BigDecimal(a);
        bigDecimal = bigDecimal.setScale(b, BigDecimal.ROUND_HALF_UP);
        return bigDecimal.doubleValue();
    }
}


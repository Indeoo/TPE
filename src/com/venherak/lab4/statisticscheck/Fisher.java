package com.venherak.lab4.statisticscheck;

import com.venherak.lab4.Matrix;

public class Fisher {
    private double[][] my;
    private double[][] mx;
    private int[] cb;
    private double[] b;
    private int d;
    private double Fp;
    private double[][] table;
    public Fisher(double[][] my, double[][] mx, int[] cb, double[] b) {
        this.my = my;
        this.cb = cb;
        this.mx = mx;
        this.b = b;
        table = new double[][]{
                {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12}, //f4
                {1, 161.45, 199.50, 215.71, 224.58, 230.16, 233.99, 236.77, 238.88, 240.54, 241.88, 245.95},
                {2, 18.51, 19.00, 19.16, 19.25, 19.30, 19.33, 19.35, 19.37, 19.38, 19.40, 19.43},
                {3, 10.13, 9.55, 9.28, 9.12, 9.01, 8.94, 8.89, 8.85, 8.81, 8.79, 8.70},
                {4, 7.71, 6.94, 6.59, 6.39, 6.26, 6.16, 6.09, 6.04, 6.00, 5.96, 5.86},
                {5, 6.61, 5.79, 5.41, 5.19, 5.05, 4.95, 4.88, 4.82, 4.77, 4.74, 4.62},
                {6, 5.99, 5.14, 4.76, 4.53, 4.39, 4.28, 4.21, 4.15, 4.10, 4.06, 3.94},
                {7, 5.59, 4.74, 4.35, 4.12, 3.97, 3.87, 3.79, 3.73, 3.68, 3.64, 3.51},
                {8, 5.32, 4.46, 4.07, 3.84, 3.69, 3.58, 3.50, 3.44, 3.39, 3.35, 3.22},
                {9, 5.12, 4.26, 3.86, 3.63, 3.48, 3.37, 3.29, 3.23, 3.18, 3.14, 3.01},
                {10, 4.96, 4.10, 3.71, 3.48, 3.33, 3.22, 3.14, 3.07, 3.02, 2.98, 2.85},
                {11, 4.84, 3.98, 3.59, 3.36, 3.20, 3.09, 3.01, 2.95, 2.90, 2.85, 2.72},
                {12, 4.75, 3.89, 3.49, 3.26, 3.11, 3.00, 2.91, 2.85, 2.80, 2.75, 2.62},
                {13, 4.67, 3.81, 3.41, 3.18, 3.03, 2.92, 2.83, 2.77, 2.71, 2.67, 2.53},
                {14, 4.60, 3.74, 3.34, 3.11, 2.96, 2.85, 2.76, 2.70, 2.65, 2.60, 2.46},
                {15, 4.54, 3.68, 3.29, 3.06, 2.90, 2.79, 2.71, 2.64, 2.59, 2.54, 2.40},
                {16, 4.49, 3.63, 3.24, 3.01, 2.85, 2.74, 2.66, 2.59, 2.54, 2.49, 2.35},
                {17, 4.45, 3.59, 3.20, 2.96, 2.81, 2.70, 2.61, 2.55, 2.49, 2.45, 2.31},
                {18, 4.41, 3.55, 3.16, 2.93, 2.77, 2.66, 2.58, 2.51, 2.46, 2.41, 2.27},
                {19, 4.38, 3.52, 3.13, 2.90, 2.74, 2.63, 2.54, 2.48, 2.42, 2.38, 2.23},
                {20, 4.35, 3.49, 3.10, 2.87, 2.71, 2.60, 2.51, 2.45, 2.39, 2.35, 2.20},
                {30, 4.2, 3.3, 2.9, 2.7, 2.5, 2.4, 2.4, 2.3, 2.3, 2.2, 2.2, 2.1},
                {120, 3.9, 3.1, 2.7, 2.5, 2.3, 2.2, 2.1, 2.05, 2.0, 1.9, 1.85, 1.8}};

    }
    public void count(){
        System.out.println();
        System.out.println("Fisher criterion");
        d = 0;
        for (int aCb : cb) {
            if (aCb == 1) {
                d++;
            }
        }
        if ( cb.length == d) {
            System.out.println("All coefficient important");
            return;
        }
        double[] disp = new double[my.length];
        double dsum = 0;
        for (int i = 0; i < disp.length; i++) {
            disp[i] = Matrix.Dispersion(my[i]);
            dsum += disp[i];
        }
        double daver = dsum / disp.length;
        double Dad = my[0].length / (my.length - d);
        dsum = 0;
        double[] yx = new double[my.length];
        for (int i = 0; i < yx.length; i++) {
            yx[i] = 0;
            for (int j = 0; j < b.length; j++) {
                if (cb[j] != 0) {
                    yx[i] += mx[i][j] * b[j];
                }
            }
        }
        for (int i = 0; i < my.length; i++) {
            dsum += Math.pow(yx[i] - Matrix.MathWaiting(my[i]), 2);
        }
        Dad *= dsum;
        System.out.println("S^2ad = " + Dad);
        Fp = Dad / daver;
        System.out.println("Fp = " + Fp);
    }
    public boolean check() {
        int f1 = my[0].length - 1;
        int f2 = my.length;
        int f3 = f2 * f1;
        int f4 = my.length - d;
        int pf3 = 0;
        int pf4 = 0;
        for (int i = 1; i < table.length; i++) {
            if (table[0][i] >= f4) {
                pf4 = i;
                break;
            }
        }
        for (int i = 1; i < table.length; i++) {
            if (table[i][0] >= f3) {
                pf3 = i;
                break;
            }
        }
        double TFp = table[pf3][pf4];
        System.out.println("TFp = " + TFp);
        if (Fp < TFp) {
            System.out.println("Fp < TFp");
            return true;
        }
        System.out.println("Fp > TFp");
        return false;
    }
}

package com.venherak.lab4;

public class Matrix {
    private double my[][];
    private double[] averY;
    private double mx[][];
    private int m;
    private double yMin;
    private double yMax;
    private double[] bcoeff;

    public Matrix(double[][] mx, double yMin, double yMax, int m) {
        this.mx = mx;
        this.yMin = yMin;
        this.yMax = yMax;
        this.m = m;

        my = new double[mx[0].length][m];
        for (int i = 0; i < my.length; i++) {
            for (int j = 0; j < m; j++) {
                my[i][j] = (yMin + Math.random() * (yMax - yMin + 1.0));
            }
        }
    }
    public double[][] getMY() {
        return my;
    }

    public void countCoeff() {
        System.out.println();
        System.out.println("Computing coefficient of equation:");

        double[] mxi = new double[mx[0].length];
        double[] mai = new double[mx[0].length];
        averY = new double[my.length];
        double y = 0;
        for (int i = 0; i < my.length; i++) {
            averY[i] = 0;
            for (int j = 0; j < my[i].length; j++) {
                averY[i] += my[i][j];
            }
            averY[i] /= my[i].length;
            y += averY[i];
        }
        for (int i = 0; i < mxi.length; i++) {
            mxi[i] = 0;
            mai[i] = 0;
            for (int j = 0; j < mx[0].length; j++) {
                mxi[i] += mx[j][i];
                mai[i] += mx[j][i] * averY[j];
            }
            mxi[i] /= mx[0].length;
            mai[i] /= mx[0].length;
        }
        for (int i = 0; i < mai.length; i++) {
            System.out.println("a" + (i + 1) + " = " + mai[i]);
        }
        double[][] maij = new double[mx[0].length][mx[0].length];
        for (int i = 0; i < maij[0].length; i++) {
            for (int j = 0; j < maij[0].length; j++) {
                maij[i][j] = 0;
                for (int g = 0; g < mx[0].length; g++) {
                    maij[i][j] += mx[g][i] * mx[g][j];
                }
                maij[i][j] /= mx[0].length;
            }
        }
        for (int i = 0; i < maij.length; i++) {
            for (int j = 0; j < maij[0].length; j++) {
                System.out.println("a" + (i + 1) + "" + (j + 1) + " = "
                        + maij[i][j]);
            }
        }
        double[][] base = new double[mx[0].length][mx[0].length];
        for (int i = 0; i < base.length; i++) {
            base[i][0] = mxi[i];
            base[0][i] = mxi[i];
        }

        for (int i = 1; i < base.length; i++) {
            for (int j = 1; j < base[i].length; j++) {
                base[i][j] = maij[i][j];
            }
        }
        double dbase = getDeterminant(base);
        double[] bi = new double[mx[0].length];
        bcoeff = new double[bi.length];
        double[][] a;
        for (int i = 0; i < bi.length; i++) {
            a = new double[base.length][base[0].length];
            for (int i2 = 0; i2 < a.length; i2++) {
                for (int j2 = 0; j2 < a[0].length; j2++) {
                    a[i2][j2] = base[i2][j2];
                }
            }
            for (int j = 0; j < mx[0].length; j++) {
                a[j][i] = mai[j];
            }
            bi[i] = getDeterminant(a);
            bcoeff[i] = bi[i] / dbase;
        }
    }

    public double[] getCoeff() {
        if (bcoeff == null) {
            countCoeff();
        }
        return bcoeff;
    }

    public double[] getAverY() {
        if (averY == null) {
            countCoeff();
        }
        return averY;
    }

    private double getDeterminant(double[][] arr) {
        double det = 0;
        int len = arr.length;
        if (len == 1) {
            det = arr[0][0];
        } else if (len == 2) {
            det = arr[0][0] * arr[1][1] - arr[0][1] * arr[1][0];
        } else
            for (int g = 0; g < len; g++) {
                double[][] minor = new double[len - 1][len - 1];
                for (int i = 0, ii = 0; i < len; i++) {
                    for (int j = 0, jj = 0; j < len; j++) {
                        if (i != 0 && g != j) {
                            minor[ii][jj] = arr[i][j];
                            jj++;
                        }
                    }
                    if (i != 0)
                        ii++;
                }
                if (g % 2 == 0) {
                    det += arr[0][g] * getDeterminant(minor);
                } else {
                    det += -arr[0][g] * getDeterminant(minor);
                }
            }
        return det;
    }

    public static double getDispersion(double[] arr) {
        double G = 0;
        double M = 0;
        for (double anArr1 : arr) {
            M += anArr1;
        }
        M = M / arr.length;
        for (double anArr : arr) {
            G += (anArr - M) * (anArr - M);
        }
        G = G / arr.length;
        return G;
    }

    public static double getExpectedValue(double[] arr) {
        double M = 0;
        for (double anArr : arr) {
            M += anArr;
        }
        M = M / arr.length;
        return M;
    }
}


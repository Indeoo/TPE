package com.venherak.lab2;

public class Experiment {

    private double R_critical = 2.10; //number = 6, p = 0.95

    private double X11 = -1.0;
    private double X12 = -1.0;
    private double X13 = +1.0;
    private double X21 = -1.0;
    private double X22 = +1.0;
    private double X23 = -1.0;

    private int minX1;
    private int maxX1;
    private int minX2;
    private int maxX2;
    private int minY;
    private int maxY;

    private double[] averageY;

    private double b0;
    private double b1;
    private double b2;

    private int numberOfY;
    private double[][] Y;

    public Experiment(int minX1, int maxX1, int minX2, int maxX2, int minY, int maxY, int numberOfY) {
        this.minX1 = minX1;
        this.maxX1 = maxX1;
        this.minX2 = minX2;
        this.maxX2 = maxX2;
        this.minY = minY;
        this.maxY = maxY;
        this.numberOfY = numberOfY;
        randomY();
        romanovskiy();
        normalize();
        naturalize();
    }

    private void randomY() {
        Y = new double[3][numberOfY];
        for (int i = 0; i < Y.length; i++) {
            for (int j = 0; j < Y[i].length; j++) {
                Y[i][j] = Math.random() * (maxY - minY) + minY;
            }
        }
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < numberOfY; j++) {
                System.out.print(Y[i][j] + "    ");
            }
            System.out.println();
        }
        System.out.println();
    }

    private void romanovskiy() {

        averageY = new double[3];
        for (int i = 0; i < numberOfY; i++) {
            averageY[0] += Y[0][i];
            averageY[1] += Y[1][i];
            averageY[2] += Y[2][i];
        }
        averageY[0] /= numberOfY;
        averageY[1] /= numberOfY;
        averageY[2] /= numberOfY;
        System.out.println("Середнє у1 = " + averageY[0]);
        System.out.println("Середнє у2 = " + averageY[1]);
        System.out.println("Середнє у3 = " + averageY[2]);
        System.out.println();

        double sigmaY1 = 0;
        double sigmaY2 = 0;
        double sigmaY3 = 0;
        for (int i = 0; i < numberOfY; i++) {
            sigmaY1 += Math.pow(Y[0][i] - averageY[0], 2);
            sigmaY2 += Math.pow(Y[1][i] - averageY[1], 2);
            sigmaY3 += Math.pow(Y[2][i] - averageY[2], 2);
        }
        sigmaY1 /= numberOfY;
        sigmaY2 /= numberOfY;
        sigmaY3 /= numberOfY;
        System.out.println("Сигма у1 = " + sigmaY1);
        System.out.println("Сигма у2 = " + sigmaY2);
        System.out.println("Сигма у3 = " + sigmaY3);
        System.out.println();

        double teta0 = Math.pow(2.0 * (2.0 * (double) numberOfY - 2.0) / (double) numberOfY / ((double) numberOfY - 4.0), 0.5);
        System.out.println("Основне відхилення тета0 = " + teta0);
        System.out.println();

        double Fuv1 = Math.max(sigmaY1, sigmaY2) / Math.min(sigmaY1, sigmaY2);
        double Fuv2 = Math.max(sigmaY1, sigmaY3) / Math.min(sigmaY1, sigmaY3);
        double Fuv3 = Math.max(sigmaY3, sigmaY2) / Math.min(sigmaY3, sigmaY2);
        System.out.println("Fuv1 = " + Fuv1);
        System.out.println("Fuv2 = " + Fuv2);
        System.out.println("Fuv3 = " + Fuv3);
        System.out.println();

        double tetauv1 = (numberOfY - 2.0) / (double) numberOfY * Fuv1;
        double tetauv2 = (numberOfY - 2.0) / (double) numberOfY * Fuv2;
        double tetauv3 = (numberOfY - 2.0) / (double) numberOfY * Fuv3;
        System.out.println("Тета uv1 = " + tetauv1);
        System.out.println("Тета uv2 = " + tetauv2);
        System.out.println("Тета uv3 = " + tetauv3);
        System.out.println();

        double Ruv1 = Math.abs(tetauv1 - 1) / teta0;
        double Ruv2 = Math.abs(tetauv2 - 1) / teta0;
        double Ruv3 = Math.abs(tetauv3 - 1) / teta0;
        System.out.println("Ruv1 = " + Ruv1);
        System.out.println("Ruv2 = " + Ruv2);
        System.out.println("Ruv3 = " + Ruv3);
        System.out.println();

        if (Ruv1 < R_critical && Ruv2 < R_critical && Ruv3 < R_critical) {
            System.out.println("Ruv < Rkp, дисперсія однорідна!");
        } else {
            System.out.println("Ruv > Rkp, дисперсія неоднорідна!");
        }
        System.out.println();
    }

    private void normalize() {
        double mx1 = (X11 + X12 + X13) / 3;
        double mx2 = (X21 + X22 + X23) / 3;
        double my = (averageY[0] + averageY[1] + averageY[2]) / 3;
        double a1 = (X11 * X11 + X12 * X12 + X13 * X13) / 3;
        double a2 = (X11 * X21 + X12 * X22 + X13 * X23) / 3;
        double a3 = (X21 * X21 + X22 * X22 + X23 * X23) / 3;
        double a11 = (X11 * averageY[0] + X12 * averageY[1] + X13 * averageY[2]) / 3;
        double a22 = (X21 * averageY[0] + X22 * averageY[1] + X23 * averageY[2]) / 3;

        double determinant = (a1 * a3 + mx1 * a2 * mx2 + mx2 * mx1 * a2 - mx2 * mx2 * a1 - mx1 * mx1 * a3 - a2 * a2);
        b0 = (my * a1 * a3 + mx1 * a2 * a22 + mx2 * a11 * a2 - a22 * a1 * mx2 - mx1 * a11 * a3 - my * a2 * a2) / determinant;
        b1 = (a3 * a11 + a22 * mx1 * mx2 + a2 * my * mx2 - mx2 * mx2 * a11 - a22 * a2 - mx1 * my * a3) / determinant;
        b2 = (a1 * a22 + a2 * mx1 * my + mx1 * mx2 * a11 - mx2 * my * a1 - mx1 * mx1 * a22 - a2 * a11) / determinant;
        System.out.println("b0 = " + b0);
        System.out.println("b1 = " + b1);
        System.out.println("b2 = " + b2);
        System.out.println();

        System.out.println("Нормоване рівняння регресії: ");
        System.out.println("y = (" + b0 + ") + (" + b1 + ") * X1 + (" + b2 + ") * X2");
        System.out.println();

        System.out.println("Перевірка:");
        System.out.println((b0 + X11 * b1 + X21 * b2) + " = " + averageY[0]);
        System.out.println((b0 + X12 * b1 + X22 * b2) + " = " + averageY[1]);
        System.out.println((b0 + X13 * b1 + X23 * b2) + " = " + averageY[2]);
        System.out.println();
    }

    private void naturalize() {
        double deltaX1 = (maxX1 - minX1) / 2;
        double deltaX2 = (maxX2 - minX2) / 2;
        double X10 = (maxX1 + minX1) / 2;
        double X20 = (maxX2 + minX2) / 2;
        double a0 = b0 - b1 * X10 / deltaX1 - b2 * X20 / deltaX2;
        double a1 = b1 / deltaX1;
        double a2 = b2 / deltaX2;

        System.out.println("a0 = " + a0);
        System.out.println("a1 = " + a1);
        System.out.println("a2 = " + a2);
        System.out.println();

        System.out.println("Натуралізоване рівняння регресії:");
        System.out.println("y = (" + a0 + ") + (" + a1 + ") * X1 + (" + a2 + ") * X2");
        System.out.println();

        System.out.println("Перевірка:");
        System.out.println((a0 + minX1 * a1 + minX2 * a2) + " = " + averageY[0]);
        System.out.println((a0 + minX1 * a1 + maxX2 * a2) + " = " + averageY[1]);
        System.out.println((a0 + maxX1 * a1 + minX2 * a2) + " = " + averageY[2]);
    }
}

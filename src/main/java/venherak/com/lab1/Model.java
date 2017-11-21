package venherak.com.lab1;

import java.util.Random;

public class Model {

    private double leftLimit;
    private double rightLimit;

    private double[] factorArray;
    private double[][] matrixPlan;
    private double[] standard;
    private double[] interval;
    private int optimum;

    Model(double a, double b, double a0, double a1, double a2, double a3) {
        leftLimit = a;
        rightLimit = b;
        factorArray = new double[4];
        factorArray[0] = a0;
        factorArray[1] = a1;
        factorArray[2] = a2;
        factorArray[3] = a3;
        this.generate();


    }

    public double[] getRow(int i) {
        return new double[]{matrixPlan[i][0], matrixPlan[i][1], matrixPlan[i][2]};
    }

    public double[][] getMatrixPlan() {
        return matrixPlan;
    }

    public double[] getStandard() {
        return standard;
    }

    public double[] getInterval() {
        return interval;
    }

    public int getOptimum() {
        return optimum;
    }

    void generate() {
        Random rand = new Random();
        matrixPlan = new double[8][7];
        for (int i = 0; i < 8; i++) {
            matrixPlan[i][3] = factorArray[0];
            for (int j = 0; j < 3; j++) {
                matrixPlan[i][j] = leftLimit + rand.nextDouble() * (rightLimit - leftLimit);
                matrixPlan[i][3] += factorArray[j + 1] * matrixPlan[i][j];
            }
        }

        double[] leftBorder = new double[3];
        double[] rightBorder = new double[3];

        for (int j = 0; j < 3; j++) {
            leftBorder[j] = matrixPlan[0][j];
            rightBorder[j] = matrixPlan[0][j];
        }

        for (int i = 1; i < 8; i++) {
            for (int j = 0; j < 3; j++) {
                if (leftBorder[j] > matrixPlan[i][j]) leftBorder[j] = matrixPlan[i][j];
                if (rightBorder[j] < matrixPlan[i][j]) rightBorder[j] = matrixPlan[i][j];
            }
        }

        standard = new double[4];
        interval = new double[3];
        for (int j = 0; j < 3; j++) {
            standard[j] = (rightBorder[j] + leftBorder[j]) / 2;
            interval[j] = (rightBorder[j] - leftBorder[j]) / 2;
            for (int i = 0; i < 8; i++) {
                matrixPlan[i][j + 4] = (matrixPlan[i][j] - standard[j]) / interval[j];
            }
        }
        standard[3] = factorArray[0];
        for (int j = 0; j < 3; j++) {
            standard[3] += factorArray[j + 1] * standard[j];
        }
        optimum = 0;
        for (int i = 1; i < 8; i++) {
            if (matrixPlan[i][3] > matrixPlan[optimum][3]) optimum = i;
        }


    }
}


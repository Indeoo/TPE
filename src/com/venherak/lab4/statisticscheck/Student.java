package com.venherak.lab4.statisticscheck;

import java.math.BigDecimal;

import static com.venherak.lab4.Matrix.getDispersion;
import static com.venherak.lab4.Matrix.getExpectedValue;
import static java.lang.Math.sqrt;

public class Student {
    private double[] t;
    private double[] table;
    private int[] coefficient;

    public Student(double[][] my, double[][] mx) {
        table = new double[] { 12.71, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365,
                2.306, 2.262, 2.228, 2.201, 2.179, 2.160, 2.145, 2.131, 2.120,
                2.110, 2.101, 2.093, 2.086, 2.080, 2.074, 2.069, 2.064, 2.060,
                2.056, 2.052, 2.048, 2.045, 2.042 };
        System.out.println("\nSTUDENT CRITERIA");
        double[] disp = new double[my.length];
        double dsum = 0;
        for (int i = 0; i < disp.length; i++) {
            disp[i] = getDispersion(my[i]);
            dsum += disp[i];
        }
        double daver = dsum / disp.length;
        System.out.println("S^2b = " + round(daver, 2));
        double dB = sqrt(daver / my.length / my[0].length);
        System.out.println("S{B} = " + round(dB, 2));
        double[] b = new double[mx[0].length];
        for (int i = 0; i < b.length; i++) {
            b[i] = 0;
            for (int j = 0; j < mx[0].length; j++) {
                b[i] += getExpectedValue(my[j]) * mx[j][i];
            }
            b[i] /= mx[0].length;
            System.out.println("B[" + i + "] = " + round(b[i], 4));
        }
        t = new double[b.length];
        for (int i = 0; i < t.length; i++) {
            t[i] = Math.abs(b[i]) / dB;
            System.out.println("t[" + i + "] = " + round(t[i], 4));
        }
        int f1 = my[0].length - 1;
        int f2 = my.length;
        int f3 = f1 * f2;
        double Tt;
        if (f3 > 30) {
            Tt = 1.960;
        } else {
            Tt = table[f3 - 1];
        }
        System.out.println("Tt = " + Tt);
        coefficient = new int[mx[0].length];
        for (int i = 0; i < coefficient.length; i++) {
            if (t[i] < Tt) {
                System.out.println("t[" + i + "] < " + Tt + " coefficient B["
                        + i + "] isn't statistically important");
                coefficient[i] = 0;
            } else {
                System.out.println("t[" + i + "] > " + Tt + " coefficient B["
                        + i + "] statistically important");
                coefficient[i] = 1;
            }
        }

    }

    public int[] getChangedCoefficient() {
        return coefficient;
    }

    private static double round(final double a, final int b) {
        BigDecimal x = new BigDecimal(a);
        x = x.setScale(b, BigDecimal.ROUND_HALF_UP);
        return x.doubleValue();
    }
}

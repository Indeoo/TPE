package com.venherak.lab4.statisticscheck;

import static com.venherak.lab4.Matrix.getDispersion;

public class Cochran {
    private double[][] my;
    private double Gp;
    private double[][] table8; // f2=N=8
    private double[][] table4;

    public Cochran(double[][] my) {
        this.my = my;
        table8 = new double[][] {
                { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 16, 36, 144, 1000 },
                { 0.6798, 0.5157, 0.4377, 0.3910, 0.3595, 0.3362, 0.3185,
                        0.3043, 0.2926, 0.2829, 0.2462, 0.2022, 0.1616, 0.1250 } };
        table4 = new double[][] {
                { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 16, 36, 144, 1000 },
                { 0.9065, 0.7679, 0.6841, 0.6287, 0.5892, 0.5598, 0.5365,
                        0.5175, 0.5017, 0.4884, 0.4366, 0.3720, 0.3093, 0.2500 } };
    }

    private void count() {
        double[] disp = new double[my.length];
        double dSum = 0;
        System.out.println();
        System.out.println("COCHRAN CRITERIA");
        System.out.println("Computing dispersion:");
        for (int i = 0; i < disp.length; i++) {
            disp[i] = getDispersion(my[i]);
            System.out.println("getDispersion[" + i + "] = " + disp[i]);
            dSum += disp[i];
        }
        double dMax = disp[0];
        for (int i = 1; i < disp.length; i++) {
            if (dMax < disp[i]) {
                dMax = disp[i];
            }
        }
        Gp = dMax / dSum;
        System.out.println("Gp = " + Gp);
    }

    public boolean check() {
        count();
        double[][] table = table8;
        if (my.length == 4) {
            table = table4;
        }
        int f1 = my[0].length - 1;
        double TGp = 0;
        for (int i = 0; i < table[0].length; i++) {
            if (!(table[0][i] < f1)) {
                TGp = table[1][i];
                System.out.println("Gt = " + TGp);
                break;
            }
        }
        if (Gp > TGp) {
            System.out.println("Gp > Gt: getDispersion isn't homogeneous. => m++");
            return false;
        } else {
            System.out.println("Gp < Gt: getDispersion is homogeneous.");
            return true;
        }
    }
}

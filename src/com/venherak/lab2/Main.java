package com.venherak.lab2;

public class Main {

    public static void main(String[] args) {
        int minX1 = -20;
        int maxX1 = 30;
        int minX2 = 30;
        int maxX2 = 80;
        int minY = -830;
        int maxY = -730;
        int numberOfY = 6;

        Experiment experiment = new Experiment(minX1, maxX1, minX2, maxX2, minY, maxY, numberOfY);
    }
}

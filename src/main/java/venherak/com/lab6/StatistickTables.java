package venherak.com.lab6;

public class StatistickTables {
    private final static double l = 1.73;

    public final static double[][] planningMatrix = {
            {-1, -1, -1},
            {-1, 1, 1},
            {1, -1, 1},
            {1, 1, -1},
            {-1, -1, 1},
            {-1, 1, -1},
            {1, -1, -1},
            {1, 1, 1},
            {-l, 0, 0},
            {l, 0, 0},
            {0, -l, 0},
            {0, l, 0},
            {0, 0, -l},
            {0, 0, l},
            {0, 0, 0}
    };

    public final static double[][] tableFisher = {
            // 1    2        3       4      5       6   12       24
            {166.4, 199.5, 215.7, 224.6, 230.2, 234.0, 244.9, 249.0},   // 01
            {18.5, 19.2, 19.2, 19.3, 19.3, 19.3, 19.4, 19.4},           // 02
            {10.1, 8.6, 9.3, 9.1, 9.0, 8.9, 8.7, 8.6},                  // 03
            {7.7, 6.9, 6.6, 6.4, 6.3, 6.2, 5.9, 5.8},                   // 04
            {6.6, 5.8, 5.4, 5.2, 5.1, 5.0, 4.7, 4.5},                   // 05
            {6.0, 5.1, 4.8, 4.5, 4.4, 4.3, 4.0, 3.8},                   // 06
            {5.5, 4.7, 4.4, 4.1, 4.0, 3.9, 3.6, 3.4},                   // 07
            {5.3, 4.5, 4.1, 3.8, 3.7, 3.6, 3.3, 3.1},                   // 08
            {5.1, 4.3, 3.9, 3.6, 3.5, 3.4, 3.1, 2.9},                   // 09
            {5.0, 4.1, 3.7, 3.5, 3.3, 3.2, 2.9, 2.7},                   // 10
            {4.8, 4.0, 3.6, 3.4, 3.2, 3.1, 2.8, 2.6},                   // 11
            {4.8, 3.9, 3.5, 3.3, 3.1, 3.0, 2.7, 2.5},                   // 12
            {4.7, 3.8, 3.4, 3.2, 3.0, 2.9, 2.6, 2.4},                   // 13
            {4.6, 3.7, 3.3, 3.1, 3.0, 2.9, 2.5, 2.3},                   // 14
            {4.5, 3.7, 3.3, 3.1, 2.9, 2.8, 2.5, 2.3},                   // 15
            {4.5, 3.6, 3.2, 3.0, 2.9, 2.7, 2.4, 2.2},                   // 16
            {4.5, 3.6, 3.2, 3.0, 2.8, 2.7, 2.4, 2.2},                   // 17
            {4.4, 3.6, 3.2, 2.9, 2.8, 2.7, 2.3, 2.1},                   // 18
            {4.4, 3.5, 3.1, 2.9, 2.7, 2.6, 2.3, 2.1},                   // 19
            {4.4, 3.5, 3.1, 2.9, 2.7, 2.6, 2.3, 2.1},                   // 20
            {4.4, 3.5, 3.1, 2.9, 2.7, 2.6, 2.3, 2.1},                   // 21
            {4.3, 3.4, 3.1, 2.8, 2.7, 2.6, 2.2, 2.0},                   // 22
            {4.3, 3.4, 3.1, 2.8, 2.7, 2.6, 2.2, 2.0},                   // 23
            {4.3, 3.4, 3.0, 2.8, 2.6, 2.5, 2.2, 2.0},                   // 24
            {4.3, 3.4, 3.0, 2.8, 2.6, 2.5, 2.2, 2.0},                   // 25
            {4.2, 3.4, 3.0, 2.7, 2.6, 2.5, 2.2, 2.0},                   // 26
            {4.2, 3.4, 3.0, 2.7, 2.6, 2.5, 2.2, 2.0},                   // 27
            {4.2, 3.3, 3.0, 2.7, 2.6, 2.4, 2.1, 1.9},                   // 28
            {4.2, 3.3, 2.9, 2.7, 2.5, 2.4, 2.1, 1.9},                   // 29
            {4.2, 3.3, 2.9, 2.7, 2.5, 2.4, 2.1, 1.9},                   // 30
            {4.1, 3.2, 2.9, 2.6, 2.5, 2.3, 2.0, 1.8},                   // 40
            {4.0, 3.2, 2.8, 2.5, 2.4, 2.3, 1.9, 1.7},                   // 60
    };

    public final static double[][] tableKohren = {
            {0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
            {0.9985, 0.9750, 0.9392, 0.9057, 0.8772, 0.8534},
            {0.9669, 0.8709, 0.7977, 0.7457, 0.7071, 0.6771},
            {0.9065, 0.7679, 0.6841, 0.6287, 0.5891, 0.5598},
            {0.8412, 0.6838, 0.5981, 0.5440, 0.5063, 0.4783},
            {0.7808, 0.6161, 0.5321, 0.4803, 0.4447, 0.4184},
            {0.7271, 0.5612, 0.4800, 0.4307, 0.3974, 0.3726},
            {0.6798, 0.5157, 0.4377, 0.3910, 0.3595, 0.3362},
            {0.6385, 0.4775, 0.4027, 0.3584, 0.3286, 0.3067},
            {0.6020, 0.4450, 0.3733, 0.3311, 0.3029, 0.2823},
            {0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
            {0.5410, 0.3924, 0.3264, 0.2880, 0.2624, 0.2439},
            {0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
            {0.4709, 0.3346, 0.2758, 0.2419, 0.2159, 0.2034}
    };

    public final static double[] tableStudent = {
            0.000, 12.71, 4.303, 3.182, 2.776, 2.571,
            2.447, 2.365, 2.306, 2.262, 2.228, 2.201,
            2.179, 2.160, 2.145, 2.131, 2.120, 2.110,
            2.101, 2.093, 2.086, 2.080, 2.074, 2.069,
            2.064, 2.060, 2.056, 2.052, 2.048, 2.045, 2.042, 1.96
    };
}
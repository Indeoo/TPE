package venherak.com.lab5;

public class Main {
    private static final int x1min = -1;
    private static final int x1max =  4;
    private static final int x2min = -3;
    private static final int x2max =  6;
    private static final int x3min = -1;
    private static final int x3max =  9;

    public static void main(String[] args) {
        Model model = new Model(x1min, x1max, x2min, x2max, x3min, x3max);
        model.calc();
        System.out.println(model.printResult());
    }
}

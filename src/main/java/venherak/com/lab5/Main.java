package venherak.com.lab5;

public class Main {

    public static void main(String[] args) {
        Model model = new Model(0, 10, -6, 6, -7, 9);
        model.calc();
        System.out.println(model.printResult());
    }
}

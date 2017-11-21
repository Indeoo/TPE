package venherak.com.lab6;

public class Main {

    public static void main(String[] args) {
        Model model = new Model();

        model.setFunkCoeff(4.3, 8.4, 6.4, 5.4, 4.1, 0.2, 7.4, 1.0, 0.3, 5.6, 2.1);
        model.calc(-20, 15, -35, 10, 10, 20);
        model.print();
    }
}


import NumericalMethodsPart1.NewtonMethod.SNonLE;

public class Main {

    public static void main(String[] args) {

        double eps1 = 10e-9, eps2 = 10e-9;
        double[] initApp = new double[]{1, 1};
        int numIt = 10;
        try {
            SNonLE.NewtonMethod(initApp,numIt,eps1,eps2);
        } catch (ArithmeticException e) {
            e.printStackTrace();
        }
    }
}
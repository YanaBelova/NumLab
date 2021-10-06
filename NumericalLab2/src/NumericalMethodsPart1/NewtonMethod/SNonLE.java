package NumericalMethodsPart1.NewtonMethod;

import GaussMethodAndFactorization.SLAE;

import java.util.Arrays;
import java.util.Scanner;

public class SNonLE {

    public static double funk1(double x1, double x2) {
        return Math.tan(x1 * x2 + 0.2) - x1 * x1;
    }

    public static double dFink1dx1(double x1, double x2) throws ArithmeticException {
        final double v = Math.cos(x1 * x2 + 0.2) * Math.cos(x1 * x2 + 0.2);
        if (v == 0) {
            throw new ArithmeticException("Деление на 0.");
        } else {
            return x2 / v - 2 * x1;
        }
    }

    public static double dFink1dx2(double x1, double x2) throws ArithmeticException {
        final double v = Math.cos(x1 * x2 + 0.2) * Math.cos(x1 * x2 + 0.2);
        if (v ==0){
            throw new ArithmeticException("Деление на 0.");
        }
        return x1 / (v);
    }

    public static double funk2(double x1, double x2) {
        return 0.5 * x1 * x1 + 2 * x2 * x2 - 1;
    }

    public static double dFink2dx1(double x1, double x2) {
        return x1;
    }

    public static double dFink2dx2(double x1, double x2) {
        return 4 * x2;
    }

    public static void NewtonMethod(double[] x, int numIt, double eps1, double eps2) throws ArithmeticException {
        int k = 1;
        double M = 0.000001;
        double[] residualVector = new double[2];
        double[][] J = new double[2][2];
        double[] tempX = new double[x.length];
        double delta1, delta2;
        double[] dx = new double[x.length];
        int value;
        System.out.println("Выберите вариант нахождения матрицы Якоби :)");
        System.out.println("Аналитический метод - 1\nКонечно-разностный метод - 2");
        Scanner in = new Scanner(System.in);
        do {
            value = in.nextInt();
        } while (value != 1 && value != 2);
        do {
//            for (int i = 0; i < x.length; i++) {
//                residualVector[i] = -x[i];
//            }
            residualVector[0] = -funk1(x[0], x[1]);
            residualVector[1] = -funk2(x[0], x[1]);
            System.out.println("Вектор невязки: " + Arrays.toString(residualVector));
            if (value == 1) {
                Jacobian(x, J);
            } else {
                Jacobian(x, J, M);
            }
            System.out.println("Якоби");
            for (double[] doubles : J) {
                System.out.println(Arrays.toString(doubles));
            }
            SLAE.GaussMethod(J, residualVector, dx);
            for (int i = 0; i < x.length; i++) {
                tempX[i] = x[i] + dx[i];
            }

            delta1 = Math.max(Math.abs(funk1(x[0], x[1])), Math.abs(funk2(x[0], x[1])));
            delta2 = Math.max(Math.abs((tempX[0] - x[0]) / ((tempX[0] >= 1) ? tempX[0] : 1)), Math.abs((tempX[1] - x[1]) / ((tempX[1] >= 1) ? tempX[1] : 1)));
            System.arraycopy(tempX, 0, x, 0, x.length);
            System.out.println("k: " + k + "\tdelta1 = " + delta1 + "\tdelta2 = " + delta2);
            if (k >= numIt) {
                throw new ArithmeticException("IER=2");
            }
            k++;
        } while (delta1 > eps1 && delta2 > eps2);
        System.out.println(Arrays.toString(x));
    }

    public static void Jacobian(double[] x, double[][] J) {
        J[0][0] = dFink1dx1(x[0], x[1]);
        J[0][1] = dFink1dx2(x[0], x[1]);
        J[1][0] = dFink2dx1(x[0], x[1]);
        J[1][1] = dFink2dx2(x[0], x[1]);
    }

    public static void Jacobian(double[] x, double[][] J, double M) {
        J[0][0] = (funk1(x[0] + M * x[0], x[1]) - funk1(x[0], x[1])) / (M * x[0]);
        J[0][1] = (funk1(x[0], x[1] + M * x[1]) - funk1(x[0], x[1])) / (M * x[1]);
        J[1][0] = (funk2(x[0] + M * x[0], x[1]) - funk2(x[0], x[1])) / (M * x[0]);
        J[1][1] = (funk2(x[0], x[1] + M * x[1]) - funk2(x[0], x[1])) / (M * x[1]);
    }
}

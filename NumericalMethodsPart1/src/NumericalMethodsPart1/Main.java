package NumericalMethodsPart1;

import java.util.Arrays;
import java.util.Scanner;
import java.util.function.Function;

public class Main {

    //---------------------------------------------------------------------------------
    public static Double funk1(Double[] x) {
        return Math.tan(x[0] * x[1] + 0.2) - x[0] * x[0];
    }

    public static Double dFink1dx1(Double[] x) throws ArithmeticException {
        final double v = Math.cos(x[0] * x[1] + 0.2) * Math.cos(x[0] * x[1] + 0.2);
        if (v == 0) {
            throw new ArithmeticException("Деление на 0.");
        } else {
            return x[1] / v - 2 * x[0];
        }
    }

    public static Double dFink1dx2(Double[] x) throws ArithmeticException {
        final double v = Math.cos(x[0] * x[1] + 0.2) * Math.cos(x[0] * x[1] + 0.2);
        if (v == 0) {
            throw new ArithmeticException("Деление на 0.");
        }
        return x[0] / (v);
    }

    public static Double funk2(Double[] x) {
        return 0.5 * x[0] * x[0] + 2 * x[1] * x[1] - 1;
    }

    public static Double dFink2dx1(Double[] x) {
        return x[0];
    }

    public static Double dFink2dx2(Double[] x) {
        return 4 * x[1];
    }

    //---------------------------------------------------------------------------------
    public static Double eulerFunction11(uNumbers values) {
        return -values.yK[0] * values.yK[1] + ((values.tK < 1e-9) ? 0.0 : (Math.sin(values.tK) / values.tK));
    }

    public static Double eulerFunction12(uNumbers values) {
        double w = 25;
        double a = 2.5 + w / 40;
        return -values.yK[1] * values.yK[1] + (a * values.tK) / (1 + values.tK * values.tK);
    }
    //---------------------------------------------------------------------------------
    public static Double eulerFunction21(uNumbers numbers) {
        return numbers.yK1[0] - numbers.yK[0] - numbers.tauK * (-numbers.yK1[0] * numbers.yK1[1] + ((numbers.tK < 1e-9) ? 0.0 : (Math.sin(numbers.tK) / numbers.tK)));
    }

    public static Double eulerFunction22(uNumbers numbers) {
        double w = 25;
        double a = 2.5 + w / 40;
        return numbers.yK1[1] - numbers.yK[1] - numbers.tauK * (-numbers.yK1[1] * numbers.yK1[1] + (a*numbers.tK) / (1 + numbers.tK * numbers.tK));
    }
    //---------------------------------------------------------------------------------


    public static final String ANSI_RESET = "\u001B[0m";
    public static final String ANSI_GREEN = "\u001B[32m";

    public static void main(String[] args) {
        Scanner in = new Scanner(System.in);
        int value;
        System.out.print("""
                Выберите номер лабораторной работы:

                Лабораторная работа №1
                РЕШЕНИЕ СИСТЕМ ЛИНЕЙНЫХ АЛГЕБРАИЧЕСКИХ УРАВНЕНИЙ
                Введите 1

                Лабораторная работа №2:
                РЕШЕНИЕ СИСТЕМ НЕЛИНЕЙНЫХ АЛГЕБРАИЧЕСКИХ УРАВНЕНИЙ
                Введите 2

                Лабораторная работа №3:
                РЕШЕНИЕ СИСТЕМ ОБЫКНОВЕННЫХ ДИФФЕРЕНЦИАЛЬНЫХ УРАВНЕНИЙ
                Введите 3

                """);
        do {
            value = in.nextInt();
        } while (value < 1 || value > 3);
        System.out.println();
        try {
            switch (value) {
                case 1 -> {
                    Double[][] matrixA = {{2.21, 3.65, 1.69, 6.99}, {8.3, 2.62, 4.1, 1.9}, {3.92, 8.45, 7.78, 2.46}, {3.77, 7.21, 8.04, 2.28}};
                    Double[][] matrixA1 = {{2.21, 3.65, 1.69, 6.99}, {8.3, 2.62, 4.1, 1.9}, {3.92, 8.45, 7.78, 2.46}, {3.77, 7.21, 8.04, 2.28}};
                    Double[][] matrixA2 = {{2.21, 3.65, 1.69, 6.99}, {8.3, 2.62, 4.1, 1.9}, {3.92, 8.45, 7.78, 2.46}, {3.77, 7.21, 8.04, 2.28}};
                    Double[] matrixB = {-8.35, -10.65, 12.21, 15.45};
                    Double[] matrixB1 = {-8.35, -10.65, 12.21, 15.45};
                    Double[] matrixX = SLAE.GaussMethod(matrixA, matrixB);
                    System.out.println(ANSI_GREEN + "Ответ:\u001B[0m" + ANSI_RESET);
                    for (int i = 0; i < matrixX.length; i++) {
                        System.out.print("x" + ++i + " = " + matrixX[--i] + "  ");
                        System.out.println();
                    }
                    Double[] residualVector = SLAE.residualVector(matrixA1, matrixB1, matrixX);
                    System.out.print("Вектор невязки: ");
                    System.out.println(Arrays.toString(residualVector));
                    System.out.println(ANSI_GREEN + "Норма: " + ANSI_RESET + SLAE.getNorm(residualVector));
                    System.out.println();
                    SLAE.RelativeError(matrixA2, matrixX);
                    System.out.println(ANSI_GREEN + "Относительная погрешность: " + ANSI_RESET + SLAE.RelativeError(matrixA2, matrixX));
                    System.out.println();
                    double l1 = 1, l2 = 10, l3 = 15;
                    Double[][] matrixALDL = {{2 * l1 + 4 * l2, 2 * (l1 - l2), 2 * (l1 - l2)}, {2 * (l1 - l2), 2 * l1 + l2 + 3 * l3, 2 * l1 + l2 - 3 * l3}, {2 * (l1 - l2), 2 * l1 + l2 - 3 * l3, 2 * l1 + l2 + 3 * l3}};
                    Double[] matrixBLDL = {-4 * l1 - 2 * l2, -4 * l1 + l2 + 9 * l3, -4 * l1 + l2 - 9 * l3};
                    if (SLAE.isSymmetric(matrixALDL)) {
                        Double[] matrixXLDL = SLAE.FactorizationMethod(matrixALDL, matrixBLDL);
                        System.out.print(ANSI_GREEN + "Ответ: " + ANSI_RESET);
                        System.out.println(Arrays.toString(matrixXLDL));
                    } else {
                        throw new ArithmeticException("Матрица не симметричная.");
                    }
                }
                case 2 -> {
                    Double[] initApp = new Double[]{1.0, 1.0};
                    Function<Double[], Double>[] functions = new Function[initApp.length];
                    functions[0] = Main::funk1;
                    functions[1] = Main::funk2;
                    Function<Double[], Double>[][] dFunctions = new Function[initApp.length][initApp.length];
                    dFunctions[0][0] = Main::dFink1dx1;
                    dFunctions[0][1] = Main::dFink1dx2;
                    dFunctions[1][0] = Main::dFink2dx1;
                    dFunctions[1][1] = Main::dFink2dx2;
                    Double[] x = SNonLE.NewtonMethod(initApp, functions, dFunctions);
                    System.out.println(ANSI_GREEN + "Ответ: " + ANSI_RESET + Arrays.toString(x));
                }
                case 3 -> {
                    Double[] u = new Double[]{0.0, -0.412};
                    Function<uNumbers, Double>[] functions1 = new Function[u.length];
                    functions1[0] = Main::eulerFunction11;
                    functions1[1] = Main::eulerFunction12;
                    Function<uNumbers, Double>[] functions2 = new Function[u.length];
                    functions2[0] = Main::eulerFunction21;
                    functions2[1] = Main::eulerFunction22;
                    System.out.println("Выберите метод решения:\nЯвный метод Эйлера → 1\nНеявный метод Эйлера → 2");
                    do {
                        value = in.nextInt();
                    } while (value < 1 || value > 2);
                    System.out.println();
                    switch (value) {
                        case 1 -> {
                            SODE.eulerExplicit(u, functions1);
                        }
                        case 2 -> {
                            SODE.eulerImplicit(u, functions2);
                        }
                    }
                }
                default -> throw new IllegalStateException("Unexpected value: " + value);
            }
        } catch (ArithmeticException e) {
            e.printStackTrace();
        }
    }
}

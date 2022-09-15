package NumericalMethodsPart1;

import java.util.function.Function;

class uNumbers {
    Double[] yK1; //y(k+1)
    Double[] yK;  //y(k)
    Double[] yK_1;//y(k-1)
    double tK;    //t(k)
    double tauK1; //tau(k+1)
    double tauK;  //tau(k)
    double tauK_1;//tau(k-1)

    public uNumbers() {
        this.yK1 = null;
        this.yK = null;
        this.yK_1 = null;
        this.tK = 0;
        this.tauK1 = 0;
        this.tauK = 0;
        this.tauK_1 = 0;
    }

    public uNumbers(Double[] yK1, Double[] yK, Double[] yK_1, double tK, double tauK1, double tauK, double tauK_1) {
        assert false;
        if (yK1.length != 0) {
            this.yK1 = new Double[yK1.length];
            System.arraycopy(yK1, 0, this.yK1, 0, yK1.length);
        }
        if (yK.length != 0) {
            this.yK = new Double[yK.length];
            System.arraycopy(yK, 0, this.yK, 0, yK.length);
        }
        if (yK_1.length != 0) {
            this.yK_1 = new Double[yK_1.length];
            System.arraycopy(yK_1, 0, this.yK_1, 0, yK_1.length);
        }
        this.tK = tK;
        this.tauK1 = tauK1;
        this.tauK = tauK;
        this.tauK_1 = tauK_1;
    }

    public uNumbers(uNumbers values) {
        assert false;
        if (values.yK1.length != 0) {
            this.yK1 = new Double[values.yK1.length];
            System.arraycopy(values.yK1, 0, this.yK1, 0, values.yK1.length);
        }
        if (values.yK.length != 0) {
            this.yK = new Double[values.yK.length];
            System.arraycopy(values.yK, 0, this.yK, 0, values.yK.length);
        }
        if (values.yK_1.length != 0) {
            this.yK_1 = new Double[values.yK_1.length];
            System.arraycopy(values.yK_1, 0, this.yK_1, 0, values.yK_1.length);
        }
        this.tK = values.tK;
        this.tauK1 = values.tauK1;
        this.tauK = values.tauK;
        this.tauK_1 = values.tauK_1;
    }
}

public class SODE {
    private static final double M = 0.0001;
    private static final double eps = 1e-3;
    private static final double tauMax = 1e-2;
    private static final double tauMin = 1e-15;
    private static final double T = 1;

    public static void eulerExplicit(Double[] u, Function<uNumbers, Double>[] functions) //явный метод Эйлера
    {
        double tau;
        uNumbers values = new uNumbers();
        values.yK = u;
        Double[] funkTempVal = new Double[values.yK.length];
        for (int count = 1; values.tK < T; count++) {
            for (int i = 0; i < values.yK.length; ++i) {
                funkTempVal[i] = functions[i].apply(values);
            }
            tau = eps / (Math.max(Math.abs(funkTempVal[0]) + eps / tauMax, Math.abs(funkTempVal[1]) + eps / tauMax));
            for (int i = 0; i < values.yK.length; ++i) { //находим значение функции при некотором параметре t
                values.yK[i] += tau * funkTempVal[i];
            }
            values.tK += tau; //делаем шаг
            println(values.yK, values.tK, count);
        }
    }

    public static void eulerImplicit(Double[] u, Function<uNumbers, Double>[] functions) {
        uNumbers values = new uNumbers(u, u, u, 0, 0, tauMin, tauMin);
        double tk1;
        Double[] epsK = new Double[values.yK.length];
        for (int count = 1; values.tK < T; count++) {
            tk1 = values.tK + values.tauK;
            NewtonMethod(values, functions);
            for (int i = 0; i < epsK.length; ++i)  // вычисляем  epsK по формуле
                epsK[i] = -(values.tauK / (values.tauK + values.tauK_1)) * (values.yK1[i] - values.yK[i] - values.tauK * (values.yK[i] - values.yK_1[i]) / values.tauK_1);
            boolean flag = false;
            for (Double aDouble : epsK)
                if (Math.abs(aDouble) > eps && !flag) //то возвращаемся начало for
                {
                    values.tauK /= 2;
                    tk1 = values.tK;
                    System.arraycopy(values.yK, 0, values.yK1, 0, values.yK.length);
                    flag = true;
                }
            if (flag)
                continue;
            Double[] temp = new Double[values.yK.length];
//            for (int i = 0; i < u.length; i++)  // определить шаг   Tau1 по формулам 3.18  и выбирается минимальным среди всех значений  Tau1
//            {
//                if (Math.abs(epsk[i]) > eps)
//                    temp[i] = values.tauK / 2;
//                if (eps / 4 < Math.abs(epsk[i]) && std::abs (epsk[i]) <= eps)
//                temp[i] = values.tauK;
//                if (Math.abs(epsk[i]) <= eps / 4)
//                    temp[i] = 2 * values.tauK;
//            }
            for (int i = 0; i < u.length; i++) {
                temp[i] = Math.sqrt(eps / Math.abs(epsK[i])) * values.tauK;
            }
            for (int i = 0; i < u.length - 1; i++) // выбирается минимальным среди всех значений  Tau1 3.19
            {
                values.tauK1 = Math.min(temp[i], temp[i + 1]);
            }
            if (values.tauK1 > tauMax)
                values.tauK1 = tauMax;
            println(values.yK, values.tK, count);
            for (int i = 0; i < u.length; i++) {
                values.yK_1[i] = values.yK[i];
                values.yK[i] = values.yK1[i];
            }
            values.tauK_1 = values.tauK;
            values.tauK = values.tauK1;
            values.tK = tk1;
        }
    }

    public static void NewtonMethod(uNumbers values, Function<uNumbers, Double>[] functions) throws ArithmeticException {
        int k = 1;
        int numIt = 100;
        double eps1 = 1e-9;
        double eps2 = 1e-9;
        Double[] residualVector = new Double[values.yK1.length];
        Double[][] J;
        Double[] tempX = new Double[values.yK1.length];
        double delta1, delta2;
        Double[] dx;
        do {
            for (int i = 0; i < residualVector.length; i++) {
                residualVector[i] = -functions[i].apply(values);
            }
            J = Jacobean(values, functions, M);
            dx = SLAE.GaussMethod(J, residualVector);
            for (int i = 0; i < values.yK1.length; i++) {
                tempX[i] = values.yK1[i] + dx[i];
            }
            delta1 = Math.max(Math.abs(functions[0].apply(values)), Math.abs(functions[1].apply(values)));
            delta2 = Math.max(Math.abs((tempX[0] - values.yK1[0]) / ((tempX[0] >= 1) ? tempX[0] : 1)), Math.abs((tempX[1] - values.yK1[1]) / ((tempX[1] >= 1) ? tempX[1] : 1)));
            System.arraycopy(tempX, 0, values.yK1, 0, values.yK1.length);
            //System.out.println("k: " + k + "\tdelta1 = " + delta1 + "\tdelta2 = " + delta2);
            if (k >= numIt) {
                throw new ArithmeticException("IER=2");
            }
            k++;
        } while (delta1 > eps1 && delta2 > eps2);
    }

    public static Double[][] Jacobean(uNumbers values, Function<uNumbers, Double>[] functions, double M) {
        Double[][] J = new Double[values.yK1.length][values.yK1.length];
        uNumbers tempValues1 = new uNumbers(values);
        uNumbers tempValues2 = new uNumbers(values);
        tempValues1.yK1[0] += M;
        tempValues2.yK1[1] += M;
        J[0][0] = (functions[0].apply(tempValues1) - functions[0].apply(values)) / M;
        J[0][1] = (functions[0].apply(tempValues2) - functions[0].apply(values)) / M;
        J[1][0] = (functions[1].apply(tempValues1) - functions[1].apply(values)) / M;
        J[1][1] = (functions[1].apply(tempValues2) - functions[1].apply(values)) / M;
        return J;
    }

    public static void println(Double[] yK, double tK, int count) {
        System.out.println("Номер итерации: " + count);
        for (int i = 0; i < yK.length; i++) {
            System.out.println("u[" + (i + 1) + "]" + " = " + yK[i]);
        }
        System.out.println("t" + count + " = " + tK);
        System.out.println();
    }
}
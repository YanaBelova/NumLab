package com.company;

import static java.lang.Math.*;

public class Main {
    static double f (double x) {
        return sqrt(1 + pow(x, 3));            //выражение с х   //вариант 1
    }
    static double f (double x, double y) {
        return 1 / (pow(x+y, 2));              //выражение с x and y   //вар 30
    }
    static void TrapeziumIntegral(double a, double b, int n, double EPS) {
        double h = (b - a) / n;
        double integral = 0.0, integral_prev = 0.0, sum = 0.0;

        do {
            integral_prev = integral;
            integral = f(a) + f(b);
            sum = 0.0;
            for (int i = 1; i < n; i++) {
                sum += 2 * f(a + i*h);
            }
            integral += sum;
            integral *= h / 2;

            h /= 2;
            n *=2;
            System.out.println( "n= " + n+ "  Интеграл в цикле = " + integral);
        } while (abs(integral_prev - integral) >= 3 * EPS);
        System.out.println( "По формуле трапеции = " + integral+"\n\n\n");
    }

    static void SimpsonIntegral(double a, double b, double EPS) {
        int m = 1, n = 2 * m;
        double h = (b - a) / n;
        double integral = 0.0, integral_prev = 0.0, sum = 0.0;

        do {
            integral_prev = integral;
            integral = f(a) + f(b);
            sum = 0;
            for (int i = 1; i < n; i++) {
                if (i % 2 != 0)
                    sum += 4 * f(a + i*h);
                else
                    sum += 2 * f(a + i*h);
            }
            integral += sum;
            integral *= h / 3;

            h /= 2;
            n *= 2;
            System.out.println( "n= " + n+ "  Интеграл в цикле = " + integral);
        } while (abs(integral_prev - integral) >= 15 * EPS);
        System.out.println( "По формуле Симпсона = " + integral + "\n\n\n");
    }


    static void SquareSimpsonIntegral(double a, double b, double c, double d, double eps2) {
        int m = 5; int n = 2;
        double sum = 0; double sum_prev = 0;

        do {
            sum_prev = sum;
            sum = 0;
            double hx = (b - a)/(2*n);
            double hy = (d - c)/(2*m);

            double xi = a;
            double yi = c;

            double[] Xi = new double[2*n + 1];
            Xi[0] = xi;

            for (int i = 1; i <= 2*n; i++)
                Xi[i] = Xi[i - 1] + hx;

            double[] Yi = new double[2*m + 1];
            Yi[0] = yi;

            for (int j = 1; j <= 2*m; j++)
                Yi[j] = Yi[j - 1] + hy;

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    sum += f(Xi[2*i], Yi[2*j]);
                    sum += 4 * f(Xi[2*i+1], Yi[2*j]);
                    sum += f(Xi[2*i+2], Yi[2*j]);
                    sum += 4 * f(Xi[2*i], Yi[2*j+1]);
                    sum += 16 * f(Xi[2*i+1], Yi[2*j+1]);
                    sum += 4 * f(Xi[2*i+2], Yi[2*j+1]);
                    sum += f(Xi[2*i], Yi[2*j+2]);
                    sum += 4 * f(Xi[2*i+1], Yi[2*j+2]);
                    sum += f(Xi[2*i+2], Yi[2*j+2]);
                }
            }
            sum *= (hx*hy/9);
            n *= 2;
            m *= 2;
            System.out.println("Cube Simpsson: " + sum + " при n = "+n + " при m = " + m);
        } while (abs(sum - sum_prev) > eps2);
        System.out.println( "По кубатурной формуле Симпсона = " + sum);
    }

    public static void main(String[] args) {

            double a = 3.0, b = 4.0, c = 1.0, d = 2.0; // интервалы интегрирования
            int n = 1; // количество разбиений фигуры
            int m = 1000; // количество разбиений фигуры
            double eps = 1E-7;
            TrapeziumIntegral(a, b, n, eps);
            SimpsonIntegral(a, b, eps);
            SquareSimpsonIntegral(a, b, c, d, eps);
    }
}

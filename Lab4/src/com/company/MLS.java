package com.company;


import java.util.Arrays;

import static java.lang.Math.*;

public class MLS {

    static void zeroMatrix(Double[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix.length; j++) {
                matrix[i][j] = 0.0;
            }
        }
    }

    //среднеквадратическое отклонение
    static Double calculateVarianceSquared(Double[] DataX, Double[] DataY, Double[] coefficientMatrix, int m, int N) {
        Double S = 0.0, temp;
        for (int i = 0; i < N; i++) {
            temp = DataY[i];
            for (int j = 0; j < m + 1; j++)
                temp -= coefficientMatrix[j] * pow(DataX[i], j);
            S += temp * temp;
        }
        S /= (N - m - 1);
        return S;
    }

    public static void leastSquareMethod(Double[] DataX, Double[] DataY, int m, int N) {
        System.out.println(Arrays.toString(DataX));
        System.out.println(Arrays.toString(DataY));
        Double[] powerX = new Double[2 * m];
        for (int i = 0; i < 2 * m; i++)
            powerX[i] = 0.0;

        for (int k = 0; k < 2 * m; k++) {
            for (int i = 0; i < N; i++) {
                powerX[k] += pow(DataX[i], k + 1);
            }
        }
        System.out.println("POWERX (Максим сказал это какие-то суммы из методы)");
        System.out.println(Arrays.toString(powerX));

        Double[][] sumX = new Double[m + 1][m + 1];    //матрица А для Гаусса
        zeroMatrix(sumX);
        for (int i = 0; i < sumX.length; i++) {
            for (int j = 0; j < sumX.length; j++) {
                if (i == 0 & j == 0) {
                    sumX[0][0] = (double) N;
                    continue;
                }
                sumX[i][j] = powerX[i + j - 1];
            }
        }
        System.out.println();
        System.out.println("SUMX (или кривое отображение)");
        for (Double[] x : sumX) {
            System.out.println(Arrays.toString(x));
        }

        Double[] praw = new Double[m + 1];      //матрица В для Гаусса
        Arrays.fill(praw, (double) 0);
        for (int k = 0; k < praw.length; k++) {
            for (int i = 0; i < N; i++)
                praw[k] += (DataY[i] * pow(DataX[i], k));
        }
        System.out.println("PRAV (правая часть как из методы)");
        System.out.println(Arrays.toString(praw));

        Double[] coefficientMatrix;
        Double[][] copySumX = new Double[m + 1][m + 1];
        for (int i = 0; i < sumX.length; i++) {
            System.arraycopy(sumX[i], 0, copySumX[i], 0, sumX.length);
        }
        coefficientMatrix = SLAE.GaussMethod(copySumX, praw);
        System.out.println(Arrays.toString(coefficientMatrix));

        System.out.println("sigma = " + sqrt(calculateVarianceSquared(DataX, DataY, coefficientMatrix, m, N)));
// насколько близки результаты к исходным

        System.out.println("v = " + coefficientMatrix[0] + " + " + coefficientMatrix[1] + "x " + "+" +coefficientMatrix[2]);
    }
}

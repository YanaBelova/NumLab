package com.NumericalLab;

import java.util.Scanner;
import java.util.Arrays;

public class NumericalLab {

    public void haha(){
        System.out.println("Max spizdil tachku");
    }

    public static boolean isSymmetric1(float mat[][]) {
        for (int i = 0; i < mat.length; i++)
            for (int j = i + 1; j < mat.length - i; j++)
                if (mat[i][j] != mat[j][i])
                    return false;
        return true;
    }

    public static float getMax(float[] inputArray) {
        float maxValue = inputArray[0];
        for (int i = 1; i < inputArray.length; i++) {
            if (inputArray[i] > maxValue) {
                maxValue = inputArray[i];
            }
        }
        return maxValue;
    }

    public static void factor(int rows, float[][] matrixA, float[] matrixB) {
        float[] matrixX = new float[rows];
        float[][] matrix_L = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        float[][] matrixA1 = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        float[][] matrixD = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        float[] matrixY = {0, 0, 0};
        float[] matrixZ = {0, 0, 0};

        for (int j = 0; j < rows; j++) {

            for (int i = j + 1; i < rows; i++) {
                if (j == 0) {
                    matrixA1[i][j] = matrixA[i][j];
                } else {
                    float temp = 0;
                    for (int k = 0; k < j; k++) {
                        temp += matrixA1[i][k] * matrix_L[j][k];
                    }
                    matrixA1[i][j] = matrixA[i][j] - temp;
                }
            }

            if (j == 0) {
                matrixD[0][0] = matrixA[0][0];
            } else {
                float temp = 0;
                for (int k = 0; k < j; k++) {
                    temp += matrixA1[j][k] * matrix_L[j][k];
                }
                matrixD[j][j] = matrixA[j][j] - temp;
            }

            for (int i = j + 1; i < rows; i++) {
                matrix_L[i][j] = matrixA1[i][j] / matrixD[j][j];
            }
        }


        matrixY[0] = matrixB[0];
        for (int i = 1; i < rows; i++) {
            float temp = 0;
            for (int k = 0; k < i; k++) {
                temp += matrix_L[i][k] * matrixY[k];
            }
            matrixY[i] = matrixB[i] - temp;
        }

        for (int i = 0; i < rows; i++) {
            matrixZ[i] = matrixY[i] / matrixD[i][i];
        }

        matrixX[rows - 1] = matrixZ[rows - 1];
        for (int i = rows - 2; i >= 0; i--) {
            float temp = 0;
            for (int k = i + 1; k < rows; k++) {
                temp += matrix_L[k][i] * matrixX[k];
            }
            matrixX[i] = matrixZ[i] - temp;
        }
        System.out.println(Arrays.toString(matrixX));
    }


    public static void main(String args[]) {

        Scanner s = new Scanner(System.in);
        System.out.println("Input number of equations: ");
        int n = s.nextInt();
        System.out.println("Input number of unknowns: ");
        int m = s.nextInt();
        float[][] A = new float[10][10];
        float[][] A1 = new float[10][10];
        float[][] A2 = new float[10][10];
        float[] b = new float[10];
        float[] b1 = new float[10];
        float[] b2 = new float[10];
        System.out.println("Input coefficients: ");
        for (int i = 0; i < n; i++) {
            A[i] = new float[10];
            for (int j = 0; j < m; j++) {
                A2[i][j] = A1[i][j] = A[i][j] = s.nextFloat();
            }
            b2[i] = b1[i] = b[i] = s.nextFloat();
        }

        int N = n;
        for (int p = 0; p < N; p++) {

            int max = p;
            for (int i = p + 1; i < N; i++) {
                if (Math.abs(A[i][p]) > Math.abs(A[max][p])) {
                    max = i;
                }
            }
            float[] temp = A[p];
            A[p] = A[max];
            A[max] = temp;
            float t = b[p];
            b[p] = b[max];
            b[max] = t;

            if (Math.abs(A[p][p]) <= 1e-10) {
                System.out.println("There are no solutions");
            }

            for (int i = p + 1; i < N; i++) {
                float alpha = A[i][p] / A[p][p];
                b[i] -= alpha * b[p];
                for (int j = p; j < N; j++) {
                    A[i][j] -= alpha * A[p][j];
                }
            }
        }
        // Обратный проход

        float[] x = new float[N];
        for (int i = N - 1; i >= 0; i--) {
            float sum = 0.0F;
            for (int j = i + 1; j < N; j++) {
                sum += A[i][j] * x[j];
            }
            x[i] = (b[i] - sum) / A[i][i];
        }

        /* Вывод результатов */

        if (n < m) {
            System.out.print("infinitely many solutions");
        } else {
            System.out.println("Solutions:");
            for (int i = 0; i < N; i++) {
                System.out.print("X" + i + " = " + x[i] + "  ");
            }
        }
        /*вычисление вектора невязки*/
        System.out.println("\n\n\n\n\nresidual vector");
        float[] result = new float[n];
        float[] resvector = new float[n];
        for (int i = 0; i < n; i++) {
            result[i] = 0;
            for (int j = 0; j < m; j++) {
                System.out.println(A1[i][j] = A1[i][j] * x[j]);
                System.out.println(result[i] = result[i] + A1[i][j]);
            }
            System.out.println("Result  " + result[i]);
            resvector[i] = b1[i] - result[i];
            System.out.println("residual vector V" + i + " = " + resvector[i]);
        }

        /*вычисление нормы этого вашего вектора невязки*/
        System.out.println("Norma: " + getMax(resvector));

        System.out.println(isSymmetric1(A2));
        if (isSymmetric1(A2)) {
            factor(n, A2, b2);
        } else {
            System.out.println("Matrix isn't simmetric");
        }
    }
}

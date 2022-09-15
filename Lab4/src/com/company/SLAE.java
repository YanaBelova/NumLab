package com.company;

import java.util.Arrays;
import java.util.Objects;

public class SLAE {

    public static Double[] GaussMethod(Double[][] matrixA, Double[] matrixB) {
        Double[] matrixX = new Double[matrixB.length];
        forwardStroke(matrixA, matrixB);          //Прямой ход
        reverseStroke(matrixA, matrixB, matrixX); //Обратный ход
        return matrixX;
    }
    //Прямой ход
    public static void forwardStroke(Double[][] matrixA, Double[] matrixB) throws ArithmeticException {
        for (int p = 0; p < matrixA.length; p++) {
            int max = p; //номер опорного элемента
            for (int i = p + 1; i < matrixA.length; i++) {
                if (Math.abs(matrixA[i][p]) > Math.abs(matrixA[max][p])) {
                    max = i;
                }
            }
            swap(matrixA[p], matrixA[max]);
            swap(matrixB[p], matrixB[max]);
            if (Math.abs(matrixA[p][p]) == 0) {
                throw new ArithmeticException("Деление на 0.");
            }
            for (int i = p + 1; i < matrixA.length; i++) {
                Double alpha = matrixA[i][p] / matrixA[p][p];
                matrixB[i] -= alpha * matrixB[p];
                for (int j = p; j < matrixA.length; j++) {
                    matrixA[i][j] -= alpha * matrixA[p][j];
                }
            }
        }
    }
    //Обратный ход
    public static void reverseStroke(Double[][] matrixA, Double[] matrixB, Double[] matrixX) {
        for (int i = matrixX.length - 1; i >= 0; i--) {
            Double sum = 0.0;
            for (int j = i + 1; j < matrixA.length; j++) {
                sum += matrixA[i][j] * matrixX[j];
            }
            matrixX[i] = (matrixB[i] - sum) / matrixA[i][i];
        }
    }

    public static void swap(Double[] A, Double[] B) {
        Double[] tempA = A;
        A = B;
        B = tempA;
    }

    public static void swap(Double A, Double B) {
        Double tempA = A;
        A = B;
        B = tempA;
    }
    //Вектор невязки
    public static Double[] residualVector(Double[][] matrixA1, Double[] matrixB1, Double[] matrixX) throws ArithmeticException {
        Double[] multiResult = multiplyMatrix(matrixA1, matrixX);
        Double[] residualVector = new Double[matrixB1.length];
        for (int i = 0; i < matrixA1.length; i++)
            residualVector[i] = multiResult[i] - matrixB1[i];
        return residualVector;
    }

    public static Double getNorm(Double[] residualVector){
        return getAbsMax(residualVector);
    }

    //Относительная погрешность
    public static Double RelativeError(Double[][] matrixA, Double[] matrixX) {
        Double[] matrixX1 = GaussMethod(matrixA, multiplyMatrix(matrixA, matrixX));
        Double[] deltaX = new Double[matrixX.length];
        for (int i = 0; i < matrixX.length; i++) {
            deltaX[i] = matrixX1[i] - matrixX[i];
        }
        return getAbsMax(deltaX) / getAbsMax(matrixX);
    }

    //Метод факторизации LDL^T
    public static Double[] FactorizationMethod(Double[][] matrixA, Double[] matrixB) {
        System.out.println("Метод LDL^T - факторизации.");
        Double[] matrixX = new Double[matrixB.length];
        Double[][] matrix_L = new Double[matrixA.length][matrixA.length];
        Double[][] matrixA1 = new Double[matrixA.length][matrixA.length];
        Double[][] matrixD = new Double[matrixA.length][matrixA.length];
        Double[] matrixY = new Double[matrixA.length];
        Double[] matrixZ = new Double[matrixA.length];
        init(matrixA1,matrix_L,matrixD,matrixY,matrixZ);
        for (int j = 0; j < matrixA1.length; j++) {
            for (int i = j + 1; i < matrixA1.length; i++) {
                if (j == 0) {
                    matrixA1[i][j] = matrixA[i][j];
                } else {
                    Double temp = (double) 0;
                    for (int k = 0; k < j; k++) {
                        temp += matrixA1[i][k] * matrix_L[j][k];
                    }
                    matrixA1[i][j] = matrixA[i][j] - temp;
                }
            }
            if (j == 0) {
                matrixD[0][0] = matrixA[0][0];
            } else {
                Double temp = (double) 0;
                for (int k = 0; k < j; k++) {
                    temp += matrixA1[j][k] * matrix_L[j][k];
                }
                matrixD[j][j] = matrixA[j][j] - temp;
            }
            for (int i = j + 1; i < matrixA1.length; i++) {
                matrix_L[i][j] = matrixA1[i][j] / matrixD[j][j];
            }
        }
        System.out.println("Матрица L:");
        for (Double[] Doubles : matrix_L) {
            System.out.println(Arrays.toString(Doubles));
        }
        System.out.println();
        System.out.println("Матрица D:");
        for (Double[] Doubles : matrixD) {
            System.out.println(Arrays.toString(Doubles));
        }
        System.out.println();
        matrixY[0] = matrixB[0];
        for (int i = 1; i < matrixA1.length; i++) {
            Double temp = (double) 0;
            for (int k = 0; k < i; k++) {
                temp += matrix_L[i][k] * matrixY[k];
            }
            matrixY[i] = matrixB[i] - temp;
        }
        for (int i = 0; i < matrixA1.length; i++) {
            matrixZ[i] = matrixY[i] / matrixD[i][i];
        }
        matrixX[matrixA1.length - 1] = matrixZ[matrixA1.length - 1];
        for (int i = matrixA1.length - 2; i >= 0; i--) {
            Double temp = (double) 0;
            for (int k = i + 1; k < matrixA1.length; k++) {
                temp += matrix_L[k][i] * matrixX[k];
            }
            matrixX[i] = matrixZ[i] - temp;
        }
        return matrixX;
    }

    public static void init(Double[][] matrixA1, Double[][]matrix_L, Double[][] matrixD, Double[] matrixY, Double[] matrixZ){
        for (int i = 0; i < matrixA1.length; i++) {
            for (int j = 0; j < matrixA1.length; j++) {
                if (i == j) {
                    matrix_L[i][j] = 1.0;
                } else {
                    matrix_L[i][j] = (double) 0;
                }
                matrixA1[i][j] = (double) 0;
                matrixD[i][j] = (double) 0;
            }
            matrixY[i] = (double) 0;
            matrixZ[i] = (double) 0;
        }
    }

    public static Double getAbsMax(Double[] inputArray) {
        double maxValue = Math.abs(inputArray[0]);
        for (int i = 1; i < inputArray.length; i++) {
            if (Math.abs(inputArray[i]) > maxValue) {
                maxValue = Math.abs(inputArray[i]);
            }
        }
        return maxValue;
    }

    public static boolean isSymmetric(Double[][] matrix) {
        for (int i = 0; i < matrix.length; i++)
            for (int j = i + 1; j < matrix.length - i; j++)
                if (!Objects.equals(matrix[i][j], matrix[j][i]))
                    return false;
        return true;
    }

    public static Double[] multiplyMatrix(Double[][] matrixA, Double[] matrixX) {
        Double[] matrixB = new Double[matrixX.length];
        Arrays.fill(matrixB, (double) 0);
        if (matrixA.length != matrixX.length) {
            throw new ArithmeticException("Матрицы не согласованы.");
        } else {
            for (int j = 0; j < matrixA.length; ++j) {
                for (int k = 0; k < matrixX.length; ++k) {
                    matrixB[j] += matrixA[j][k] * matrixX[k];
                }
            }
        }
        return matrixB;
    }
}
package GaussMethodAndFactorization;

import java.util.Arrays;

public class SLAE {
    public static final String ANSI_RESET = "\u001B[0m";
    public static final String ANSI_GREEN = "\u001B[32m";

    public static void GaussMethod(double[][] matrixA, double[] matrixB, double[] matrixX) throws ArithmeticException {
//Прямой ход
        System.out.println("Метод Гаусса.");
        System.out.println("Треугольная матрица:");
        for (int p = 0; p < matrixA.length; p++) {
            int max = p; //номер опорного элемента
            for (int i = p + 1; i < matrixA.length; i++) {
                if (Math.abs(matrixA[i][p]) > Math.abs(matrixA[max][p])) {
                    max = i;
                }
            }
            double[] temp = matrixA[p];
            matrixA[p] = matrixA[max];
            matrixA[max] = temp;
            double t = matrixB[p];
            matrixB[p] = matrixB[max];
            matrixB[max] = t;
            if (Math.abs(matrixA[p][p]) == 0) {
                throw new ArithmeticException("Деление на 0.");
            }
            for (int i = p + 1; i < matrixA.length; i++) {
                double alpha = matrixA[i][p] / matrixA[p][p];
                matrixB[i] -= alpha * matrixB[p];
                for (int j = p; j < matrixA.length; j++) {
                    matrixA[i][j] -= alpha * matrixA[p][j];
                }
            }
            System.out.println(Arrays.toString(matrixA[p]));
        }
        System.out.println();
// Обратный ход
        for (int i = matrixX.length - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < matrixA.length; j++) {
                sum += matrixA[i][j] * matrixX[j];
            }
            matrixX[i] = (matrixB[i] - sum) / matrixA[i][i];
        }
//Вывод результатов
        System.out.println(ANSI_GREEN + "Ответ:\u001B[0m" + ANSI_RESET);
        for (int i = 0; i < matrixX.length; i++) {
            System.out.print("x" + ++i + " = " + matrixX[--i] + " ");
            System.out.println();
        }
        System.out.println();
    }

    public static void residualVector(double[][] matrixA1, double[] matrixB1, double[] matrixX, double[] residualVector) throws ArithmeticException {
        System.out.print("Вектор невязки: ");
        double[] multiResult = new double[matrixX.length];
        for (int i = 0; i < matrixA1.length; i++) {
            multiResult[i] = 0;
            for (int j = 0; j < matrixA1.length; j++) {
                matrixA1[i][j] = matrixA1[i][j] * matrixX[j];
                multiResult[i] = multiResult[i] + matrixA1[i][j];
            }
            residualVector[i] = multiResult[i] - matrixB1[i];
        }
        System.out.println(Arrays.toString(residualVector));
        double norm = getAbsMax(residualVector);
        System.out.println(ANSI_GREEN + "Норма: " + ANSI_RESET + norm);
        System.out.println();
//return norm;
    }

    //Относительная погрешность
    public static void RelativeError(double[][] matrixA, double[] matrixX) {
        double[] matrixX1 = new double[matrixX.length];
        GaussMethod(matrixA, multiplyMatrix(matrixA, matrixX), matrixX1);
        double[] deltaX = new double[matrixX.length];
        for (int i = 0; i < matrixX.length; i++) {
            deltaX[i] = matrixX1[i] - matrixX[i];
        }
        double delta = getAbsMax(deltaX) / getAbsMax(matrixX);
        System.out.println(ANSI_GREEN + "Относительная погрешность: " + ANSI_RESET + delta);
        System.out.println();
    }

    //Метод факторизации LDL^T
    public static void FactorizationMethod(double[][] matrixA, double[] matrixB) {
        System.out.println("Метод LDL^T - факторизации.");
        double[] matrixX = new double[matrixB.length];
        double[][] matrix_L = new double[matrixA.length][matrixA.length];
        double[][] matrixA1 = new double[matrixA.length][matrixA.length];
        double[][] matrixD = new double[matrixA.length][matrixA.length];
        double[] matrixY = new double[matrixA.length];
        double[] matrixZ = new double[matrixA.length];
        for (int i = 0; i < matrixA1.length; i++) {
            for (int j = 0; j < matrixA1.length; j++) {
                if (i == j) {
                    matrix_L[i][j] = 1;
                } else {
                    matrix_L[i][j] = 0;
                }
                matrixA1[i][j] = 0;
                matrixD[i][j] = 0;
            }
            matrixY[i] = 0;
            matrixZ[i] = 0;
        }
        for (int j = 0; j < matrixA1.length; j++) {
            for (int i = j + 1; i < matrixA1.length; i++) {
                if (j == 0) {
                    matrixA1[i][j] = matrixA[i][j];
                } else {
                    double temp = 0;
                    for (int k = 0; k < j; k++) {
                        temp += matrixA1[i][k] * matrix_L[j][k];
                    }
                    matrixA1[i][j] = matrixA[i][j] - temp;
                }
            }
            if (j == 0) {
                matrixD[0][0] = matrixA[0][0];
            } else {
                double
                        temp = 0;
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
        for (double[] doubles : matrix_L) {
            System.out.println(Arrays.toString(doubles));
        }
        System.out.println();
        System.out.println("Матрица D:");
        for (double[] doubles : matrixD) {
            System.out.println(Arrays.toString(doubles));
        }
        System.out.println();
        matrixY[0] = matrixB[0];
        System.out.print(ANSI_GREEN + "Ответ: " + ANSI_RESET);
        for (int i = 1; i < matrixA1.length; i++) {
            double temp = 0;
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
            double temp = 0;
            for (int k = i + 1; k < matrixA1.length; k++) {
                temp += matrix_L[k][i] * matrixX[k];
            }
            matrixX[i] = matrixZ[i] - temp;
        }
        System.out.println(Arrays.toString(matrixX));
    }

    public static double getAbsMax(double[] inputArray) {
        double maxValue = Math.abs(inputArray[0]);
        for (int i = 1; i < inputArray.length; i++) {
            if (Math.abs(inputArray[i]) > maxValue) {
                maxValue = Math.abs(inputArray[i]);
            }
        }
        return maxValue;
    }

    public static boolean isSymmetric(double[][] matrix) {
        for (int i = 0; i < matrix.length; i++)
            for (int j = i + 1; j < matrix.length - i; j++)
                if (matrix[i][j] != matrix[j][i])
                    return false;
        return true;
    }

    public static double[] multiplyMatrix(double[][] matrixA, double[] matrixX) {
        if (matrixA.length != matrixX.length) {
            throw new ArithmeticException("Матрицы не согласованы.");
        } else {
            double[] matrixB = new double[matrixX.length];
            for (int j = 0; j < matrixA.length; ++j) {
                for (int k = 0; k < matrixX.length; ++k) {
                    matrixB[j] += matrixA[j][k] * matrixX[k];
                }
            }
            return matrixB;
        }
    }
}


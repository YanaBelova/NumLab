package com.company;

import static java.lang.Math.*;

public class Main {
    //степень аппроксимирующего полинома; число измерений
    //var 3
    static final int m = 2; static final int N = 11;
    static void enterMatrix_X(double[] A)
    {
        A[0]=0; A[1]=1; A[2]=2; A[3]=3; A[4]=4; A[5]=5;
        A[6]=6; A[7]=7; A[8]=8; A[9]=9; A[10]=10;
    }

    static void enterMatrix_Y(double[] A)
    {
        A[0]=3; A[1]=87; A[2]=156; A[3]=210; A[4]=238; A[5]=252;
        A[6]=239; A[7]=211; A[8]=158; A[9]=90; A[10]=-5;
    }

    static double[][] initial(int n)
    {
        double A[][] = new double[n][n + 1];
//        for (int i = 0; i < n; i++) {
//            A[i] = new double[n + 1];
//        }
        return A;
    }

    static void out(double[][] A, int n)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < (n+1); j++)
                System.out.println(A[i][j]);
        }
    }

    static void copyMatrix(double[][] A, double[][] copyA, int n)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n + 1; j++)
                copyA[i][j] = A[i][j];
        }
    }

    static void zeroMatrix(double[][] A, int n)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < (n+1); j++)
            {
                A[i][j] = 0;
            }
        }
    }

    static void Answer(double[] A, int n) //view answer
    {
        for (int i = 0; i < n; i++)
            System.out.println( "a[" + i + "] = " + A[i]);
    }

    static double[] gauss(double[][] matrix, int n, int m)
    {
        //prjamoj
        double elem;
        for (int j = 0; j < n; j++)
        {
            double max = 0;
            int coord_str = 0;
            for (int t = j; t < n; t++)
            {
                if (abs(matrix[t][j]) > max)
                {
                    max = abs(matrix[t][j]); coord_str = t;
                }
            }
            if (max > abs(matrix[j][j]))
            {
                double[] ptr = matrix[j];
                matrix[j] = matrix[coord_str];
                matrix[coord_str] = ptr;
            }
            elem = matrix[j][j];
            for (int c = j; c < m; c++)
            {
                matrix[j][c] /= elem;   //delenije stroki na elem
            }

            for (int i2 = j + 1; i2 < n; i2++)
            {
                elem = matrix[i2][j];
                for (int k = j; k < m; k++)
                    matrix[i2][k] -= elem * matrix[j][k];
            }
            //		out(matrix, n, m);

        }
        //obratnyj
        double[] xx = new double[m];
        xx[n - 1] = matrix[n - 1][n];
        for (int i = n - 2; i >= 0; i--)
        {
            xx[i] = matrix[i][n];
            for (int j = i + 1; j < n; j++)
                xx[i] -= matrix[i][j] * xx[j];
        }

       // for (int i = 0; i < n; i++)
         //   System.out.println( xx[i] + " ");

        return xx;
    }


    //среднеквадратическое отклонение
    static double kvadr_dispersii(double[] Datax, double[] Datay, double[] An, double n)
    {
        double S = 0, temp;
        for (int i = 0; i < N; i++)
        {
            temp = Datay[i];
            for (int j = 0; j < n; j++)
                temp -= An[j] * pow(Datax[i], j);
            S += temp*temp;
        }
        S /= (N - m - 1);
        return S;
    }

    static void main1()
    {
        double[] Datax = new double[N];
        double[] Datay = new double[N];
//	ifstream inx("fileX.txt");
//	ifstream iny("fileY.txt");
        enterMatrix_X(Datax);
        enterMatrix_Y(Datay);

//        for (int j = 0; j < N; j++)
//            System.out.println( Datax[j] + "   " + Datay[j]);

        double[] POWERX = new double[2*m];
        for (int i = 0; i < 2*m; i++)
            POWERX[i] = 0;

        for (int k = 1; k <= 2*m; k++)
            for (int i = 0; i < N; i++)
                POWERX[k-1] += pow(Datax[i], k);

        for (int i = 0; i < 2*m; i++)
            System.out.println( POWERX[i] + " ");

        double[][] SUMX = initial(m+1);
        zeroMatrix(SUMX, m + 1);

        for (int i = 1; i <= m + 1; i++) // matrica SUMX
        {
            for (int j = 1; j <= m + 1; j++)
                if ((i+j)>=3)
                    SUMX[i-1][j-1] = POWERX[i+j-3];
        }
        SUMX[0][0] = N;

        out(SUMX, m+1);
        for (int k = 0; k <= m; k++) // правые части
        {
            for (int i = 0; i < N; i++)
                SUMX[k][m+1] += (Datay[i] * pow(Datax[i], k));
        }

        out(SUMX, m+1);
        double[] An = new double[m+1];
        double[][] copyA = initial(m+1);
        copyMatrix(SUMX, copyA, m+1);

        An = gauss(copyA, m+1, m+2);
        Answer(An, m+1);

        System.out.println( "sigma = " + sqrt(kvadr_dispersii(Datax, Datay, An, m + 1)));
        // насколько близки результаты к исходным

        System.out.println( "y = " + An[2] + "x^2 + " + An[1] + "x + " + An[0]);
    }




















    static Double[] enterMatrixX(int N) {
        Double[] matrixX = new Double[N];
        matrixX[0] = 0.0;
        matrixX[1] = 1.0;
        matrixX[2] = 2.0;
        matrixX[3] = 3.0;
        matrixX[4] = 4.0;
        matrixX[5] = 5.0;
        matrixX[6] = 6.0;
        matrixX[7] = 7.0;
        matrixX[8] = 8.0;
        matrixX[9] = 9.0;
        matrixX[10] = 10.0;
        return matrixX;
    }

    static Double[] enterMatrixY(int N) {
        Double[] matrixY = new Double[N];
        matrixY[0] = 3.0;
        matrixY[1] = 87.0;
        matrixY[2] = 156.0;
        matrixY[3] = 210.0;
        matrixY[4] = 238.0;
        matrixY[5] = 252.0;
        matrixY[6] = 239.0;
        matrixY[7] = 211.0;
        matrixY[8] = 158.0;
        matrixY[9] = 90.0;
        matrixY[10] = -5.0;
        return matrixY;
    }

    public static void main(String[] args) {
      // Main.main1();
        int N = 11;
        int degreeOfApproximation = 2;        //степень аппроксимации
        System.out.println("Исходные данные:");
        Double[] matrixX = enterMatrixX(N);
        Double[] matrixY = enterMatrixY(N);
        MLS.leastSquareMethod(matrixX,matrixY,degreeOfApproximation,N);
    }
}

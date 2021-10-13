package lab5;
import java.lang.*;
public class Lab5 {
    public static void main(String[] args) {
        System.out.println("Метод трапеции: "+integralTrapec(0.8, 1.762, x -> (Math.pow((1 + x * x * x), 1 / 2))));
        Simpson(0.8, 1.762, 0.00001);
    }

    public static double integralTrapec(double a, double b, Function function) {
        double area = 0;
        int n = 1;//
        double h = (b - a) / n;//0.01;
        for (int i = 0; i < (b - a) / h; i++) {
            area += h * (0.5 * (function.func(a + i * h) + function.func(a + (i + 1) * h)));
            n++;
        }
        return area;
    }

    public static double f (double x) {
        return (Math.pow((1 + x * x * x), 1 / 2));
    }
    public static void Simpson(double a, double b, double eps){
        //Scanner in = new Scanner(System.in);

        //double a = in.nextDouble();//Нижний и верхний пределы интегрирования (a, b), погрешность (eps).
        //double b = in.nextDouble();
        //double eps = in.nextDouble();

        double I=eps+1, I1=0;//I-предыдущее вычисленное значение интеграла, I1-новое, с большим N.

        for (int N=2; (N<=4)||(Math.abs(I1-I)>eps); N*=2){
            double h, sum2=0, sum4=0, sum=0;
            h=(b-a)/(2*N);//Шаг интегрирования.
            for (int i=1; i<=2*N-1; i+=2){
                sum4+=f(a+h*i);//Значения с нечётными индексами, которые нужно умножить на 4.
                sum2+=f(a+h*(i+1));//Значения с чётными индексами, которые нужно умножить на 2.
            }
            sum=f(a)+4*sum4+2*sum2-f(b);//Отнимаем значение f(b) так как ранее прибавили его дважды.
            I=I1;
            I1=(h/3)*sum;
        }
        System.out.println("Метод Симпсона: "+I1);
    }
}
interface Function{
    public double func (double x);
}

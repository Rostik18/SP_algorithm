using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace SP_algorithm
{
    class SP_algorithm
    {
        public int n;
        public double[,] A;
        public double[,] L;
        public double[,] U;
        public double[] y0;
        public double[] y1;
        public double[] x;
        public double lambda0 = 0;
        public double lambda1 = 0;
        public double s;
        public double norma_y;

        public const double epselon = 0.0001;

        public void ReadFile()
        {
            string path = @"C:\Users\Ростик\source\repos\SP algorithm\SP algorithm\File1.txt";
            StreamReader file = new StreamReader(path);

            n = Convert.ToInt32(file.ReadLine());
            A = new double[n, n];
            L = new double[n, n];
            U = new double[n, n];
            y1 = new double[n];
            y0 = new double[n];
            x = new double[n];
            string line;
            string[] lineArr;

            for (int i = 0; i < n; i++)
            {
                line = file.ReadLine();
                lineArr = line.Split(' ');

                for (int j = 0; j < n; j++)
                {
                    A[i, j] = Convert.ToInt32(lineArr[j]);
                }
            }
            file.Close();

            for (int i = 0; i < n; i++)
                y0[i] = y1[i] = 1;

            s = Skal_Dob(y1, y1);
            norma_y = Math.Sqrt(s);
            for (int i = 0; i < n; i++)
                x[i] = y1[i] / norma_y;

            Console.Write($"iтерацiя 0: {lambda1:f5} ; ");
            for (int i = 0; i < n; i++)
                Console.Write($"{x[i]:f5} ");
            Console.WriteLine();

            LU_decay();
        }

        void Algorithm()
        {
            double t;
            int iterationNum = 0;
            double sub = 2 * epselon;
            while (sub > epselon)
            {
                Console.Write($"iтерацiя {++iterationNum}: ");

                //y1 = Matr_na_Vect(A, x);
                //y1 = LU_Solution(x);
                y1 = Gauss_Method(A, x);

                s = Skal_Dob(y1, y1);
                t = Skal_Dob(y1, x);
                norma_y = Math.Sqrt(s);

                for (int i = 0; i < n; i++)
                {
                    x[i] = y1[i] / norma_y;
                    y0[i] = y1[i];
                }

                lambda1 = s / t;

                Console.Write($"{lambda1:f5} ; ");
                for (int i = 0; i < n; i++)
                    Console.Write($"{x[i]:f5} ");
                Console.WriteLine();

                sub = lambda1 - lambda0;
                lambda0 = lambda1;
            }

        }
        public double[] Gauss_Method(double[,] _B, double[] _b)
        {
            double[] rez = new double[n];
            double[] b = new double[n];
            double[,] _A = new double[n, n];
            for (int i = 0; i < n; i++)
            {
                b[i] = _b[i];
                for (int j = 0; j < n; j++)
                {
                    _A[i, j] = _B[i, j];
                }
            }


            for (int k = 0; k < n - 1; k++)
            {
                if (_A[k, k] == 0) SwapLines(ref _A, ref b, k);

                for (int i = k + 1; i < n; i++)
                {
                    double m = -_A[i, k] / _A[k, k];
                    b[i] = b[i] + m * b[k];
                    for (int j = k + 1; j < n; j++)
                    {
                        _A[i, j] = _A[i, j] + m * _A[k, j];
                    }
                    _A[i, k] = 0;
                }
            }
            rez[n - 1] = b[n - 1] / _A[n - 1, n - 1];
            for (int k = n - 2; k >= 0; k--)
            {
                rez[k] = (b[k] - Sum(_A, k, rez)) / _A[k, k];
            }
            return rez;


            double Sum(double[,] __A, int k, double[] _x)
            {
                double _rez = 0;
                for (int j = k + 1; j < n; j++)
                {
                    _rez += __A[k, j] * _x[j];
                }
                return _rez;
            }
            void SwapLines(ref double[,] __A, ref double[] __b, int k)
            {
                int goodLine = k;
                double max = __A[k, k];
                for (int i = k; i < n; i++)
                {
                    if (__A[i, i] != 0 && i != k && Math.Abs(__A[i, i]) > Math.Abs(max))
                    {
                        max = __A[i, i];
                        goodLine = i;
                    }
                }
                for (int j = 0; j < n; j++)
                {
                    double swap_tmp = __A[k, j];
                    __A[k, j] = __A[goodLine, j];
                    __A[goodLine, j] = swap_tmp;
                }
                double tmp_b = __b[k];
                __b[k] = __b[goodLine];
                __b[goodLine] = tmp_b;
            }
        }

        double[] Matr_na_Vect(double[,] A, double[] y)
        {
            double[] rez = new double[n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    rez[i] += A[i, j] * y[j];
                }
            }
            return rez;
        }

        double Skal_Dob(double[] _y1, double[] _y2)
        {
            double rez = 0.0;
            for (int i = 0; i < _y1.Length; i++)
            {
                rez += _y1[i] * _y2[i];
            }
            return rez;
        }


        private double[] LU_Solution(double[] _b)
        {
            double[] __y = new double[n];
            double[] __b = new double[n];
            double[] __x = new double[n];
            for (int i = 0; i < n; i++)
            {
                __b = _b;
            }
            for (int i = 0; i < n; i++)
            {
                __y[i] = (__b[i] - sumY(i)) / L[i, i];
            }
            for (int i = n - 1; i >= 0; i--)
            {
                __x[i] = (__y[i] - sumX(i)) / U[i, i];
            }

            return __x;

            double sumX(int i)
            {
                double rez = 0;
                for (int k = i + 1; k < n; k++)
                {
                    rez += U[i, k] * __x[k];
                }
                return rez;
            }
            double sumY(int i)
            {
                double rez = 0;
                for (int k = 0; k < i; k++)
                {
                    rez += L[k, i] * __y[k];
                }
                return rez;
            }
        }
        private void LU_decay()
        {
            U = A;

            for (int i = 0; i < n; i++)
                for (int j = i; j < n; j++)
                    L[j, i] = U[j, i] / U[i, i];

            for (int k = 1; k < n; k++)
            {
                for (int i = k - 1; i < n; i++)
                    for (int j = i; j < n; j++)
                        L[j, i] = U[j, i] / U[i, i];

                for (int i = k; i < n; i++)
                    for (int j = k - 1; j < n; j++)
                        U[i, j] = U[i, j] - L[i, k - 1] * U[k - 1, j];
            }
        }

        void Show()
        {
            Console.WriteLine("Матриця:");
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    Console.Write(L[i, j] + " ");
                }
                Console.WriteLine();
            }
        }

        public void Run()
        {
            ReadFile();
            Algorithm();

            //show();
        }
    }
    class Program
    {
        static void Main(string[] args)
        {
            SP_algorithm algorithm = new SP_algorithm();
            algorithm.Run();

            Console.ReadKey();
        }
    }
}

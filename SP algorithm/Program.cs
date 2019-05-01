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
        public double[] y;
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
            y = new double[n];
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
                y[i] = 1;

            s = Skal_Dob(y, y);
            norma_y = Math.Sqrt(s);
            for (int i = 0; i < n; i++)
                x[i] = y[i] / norma_y;
        }

        void Algorithm()
        {
            int iterationNum = 0;
            double sub = 2 * epselon;
            while (sub > epselon)
            {
                Console.Write($"iтерацiя {++iterationNum}: ");
                //y = Matr_na_Vect(A, x);
                y = Gauss_Method(A, x);

                s = Skal_Dob(y, y);
                double t = Skal_Dob(y, x);
                norma_y = Math.Sqrt(s);
                for (int i = 0; i < n; i++)
                    x[i] = y[i] / norma_y;
                lambda1 = s / t;

                Console.Write($"{lambda1:f4} ; ");
                for (int i = 0; i < n; i++)
                    Console.Write($"{x[i]:f5} ");
                Console.WriteLine();


                sub = lambda1 - lambda0;
                lambda0 = lambda1;

            }

        }
        public double[] Gauss_Method(double[,] _B,  double[] _b)
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
        }
        private double Sum(double[,] _A, int k, double[] _x)
        {
            double rez = 0;
            for (int j = k + 1; j < n; j++)
            {
                rez += _A[k, j] * _x[j];
            }
            return rez;
        }
        private void SwapLines(ref double[,] _A, ref double[] b, int k)
        {
            int goodLine = k;
            double max = _A[k, k];
            for (int i = k; i < n; i++)
            {
                if (_A[i, i] != 0 && i != k && Math.Abs(_A[i, i]) > Math.Abs(max))
                {
                    max = _A[i, i];
                    goodLine = i;
                }
            }
            for (int j = 0; j < n; j++)
            {
                double swap_tmp = _A[k, j];
                _A[k, j] = _A[goodLine, j];
                _A[goodLine, j] = swap_tmp;
            }
            double tmp_b = b[k];
            b[k] = b[goodLine];
            b[goodLine] = tmp_b;
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

        double Skal_Dob(double[] y1, double[] y2)
        {
            double rez = 0;
            for (int i = 0; i < y1.Length; i++)
            {
                rez += y1[i] * y2[i];
            }
            return rez;
        }


        public void Run()
        {
            ReadFile();
            Algorithm();
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

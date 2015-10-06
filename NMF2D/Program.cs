using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using InoueLab;
using CsvFileIO;
using MySignalProcessing;


namespace NMF2D
{
    class Program
    {

        static double[,] X;
        static List<double[,]> T = new List<double[,]>();
        static List<double[,]> V = new List<double[,]>();
        static double[,] Xhat;
        static double[,] ones;


        static int tau; //時間ピッチ数
        static int fai; //周波数ピッチ数
        static int K;   //基底数
        static int I;
        static int J;

        static double min;

        static void Main(string[] args)
        {
            init();
            int itteration = 0;
            while (itteration < 10)
            {
                itteration++;
                Updates();
                Console.WriteLine("itteration : " + itteration);
            }
            CsvFileIO.CsvFileIO.WriteData("out_K.csv", T[0]);
        }
        //----------------------------------------------------------------------------
        static void init()
        {
            string input_filename=@"C:\Users\優\Desktop\音素材\cq.csv";
            X = CsvFileIO.CsvFileIO.ReadData(input_filename);

            I = X.GetLength(0);
            J = X.GetLength(1);
            K = 3;
            tau = 7;
            fai = 12;
            //T.Clear();

            //T,Vの初期化
            for (int tt = 0; tt < tau; tt++)
                T.Add(initMatrix(I, K));
            for (int ff = 0; ff < fai; ff++)
                V.Add(initMatrix(K, J));

            //onesの初期化
            ones = new double[I, J];
            for (int i = 0; i < I; i++)
                for (int j = 0; j < J; j++)
                    ones[i, j] = 1;


        }
        //----------------------------------------------------------------------------
        static void Updates()
        {
            updateXhat();
            updateV();
            updateXhat();
            updateT();
        }

        //----------------------------------------------------------------------------
        static void updateT()
        {
            double[,] A = new double[I, K];
            double[,] B = new double[I, K];
            double[,] XdiviXhat = new double[I, J];

            for (int i = 0; i < I; i++)
                for (int j = 0; j < J; j++)
                    if (Xhat[i, j] != 0)
                        XdiviXhat[i, j] = X[i, j] / Xhat[i, j];

            for (int tt = 0; tt < tau; tt++)
            {
                for (int ff = 0; ff < fai; ff++)
                {
                    A = Mt.Add(A, Mt.Mul(shift_up(XdiviXhat, ff), shift_right(V[ff], tt).T()));
                    B = Mt.Add(B, Mt.Mul(ones, shift_right(V[ff], tt).T()));
                }
                for (int i = 0; i < I; i++)
                    for (int k = 0; k < K; k++)
                    {
                        if (A[i, k] == 0 && B[i, k] == 0)
                            T[tt][i, k] = T[tt][i, k];
                        else if (B[i, k] != 0)
                            T[tt][i, k] = T[tt][i, k] * A[i, k] / B[i, k];
                    }
            }
        }
        //----------------------------------------------------------------------------
        static void updateV()
        {
            double[,] A = new double[K, J];
            double[,] B = new double[K, J];
            double[,] XdiviXhat = new double[I, J];

            for (int i = 0; i < I; i++)
                for (int j = 0; j < J; j++)
                    if (Xhat[i, j] != 0)
                        XdiviXhat[i, j] = X[i, j] / Xhat[i, j];

            for (int ff = 0; ff < fai; ff++)
            {
                for (int tt = 0; tt < tau; tt++)
                {
                    A = Mt.Add(A,Mt.Mul(shift_down(T[tt], ff).T(), shift_left(XdiviXhat, tt)));
                    B = Mt.Add(B,Mt.Mul(shift_down(T[tt], ff).T(), ones));
                }
                for (int k = 0; k < K; k++)
                    for (int j = 0; j < J; j++)
                    {
                        if (A[k, j] == 0 && B[k, j] == 0)
                            V[ff][k, j] = V[ff][k, j];
                        else if (B[k, j] != 0)
                            V[ff][k, j] = V[ff][k, j] * A[k, j] / B[k, j];
                    }

            }
        }
        //--------------------------------------------------------------------------------------------------
        static void updateXhat()
        {
            double[,] A;
            double[,] Sum = new double[I, J];
            for (int tt = 0; tt < tau; tt++)
            {
                for (int ff = 0; ff < fai; ff++)
                {
                    A = Mt.Mul(shift_down(T[tt], ff), shift_right(V[ff], tt));
                    Sum = Mt.Add(Sum, A);
                }
            }
            Xhat = Sum;
        }
        //--------------------------------------------------------------------------------------------------

        #region パーツ
        //--------------------------------------------------------------------------------------------------
        //初期値をランダムに設定 (0,1]
        public static double[,] initMatrix(int N, int M)
        {
            RandomMT rand = new RandomMT();
            double[,] A = new double[N, M];
            for (int i = 0; i < N; i++)
                for (int j = 0; j < M; j++)
                    A[i, j] = 10 * rand.Double32OC();
            return A;
        }
        //----------------------------------------------------------------------------
        static double[,] shift_left(double[,] m, int shift)
        {
            double[,] matrix = new double[m.GetLength(0), m.GetLength(1)];
            if (shift < m.GetLength(1))
            {
                for (int i = 0; i < m.GetLength(0); i++)
                {
                    for (int j = 0; j < m.GetLength(1) - shift; j++)
                        matrix[i, j] = m[i, j + shift];
                    for (int j = m.GetLength(1) - shift; j < m.GetLength(1); j++)
                        matrix[i, j] = 0;
                }
            }
            else Console.WriteLine("shift spot error");
            return matrix;
        }
        //----------------------------------------------------------------------------
        static double[,] shift_right(double[,] m, int shift)
        {
            double[,] matrix = new double[m.GetLength(0), m.GetLength(1)];
            if (shift < m.GetLength(1))
            {
                for (int i = 0; i < m.GetLength(0); i++)
                {
                    for (int j = m.GetLength(1) - shift - 1; j >= 0; j--)
                        matrix[i, j + shift] = m[i, j];
                    for (int j = 0; j < shift; j++)
                        matrix[i, j] = 0;
                }
            }
            else Console.WriteLine("shift spot error");
            return matrix;
        }
        //----------------------------------------------------------------------------------------
        static double[,] shift_down(double[,] m, int shift)
        {
            double[,] matrix = new double[m.GetLength(0), m.GetLength(1)];
            if (shift < m.GetLength(0))
            {
                for (int j = 0; j < m.GetLength(1); j++)
                {
                    for (int i = m.GetLength(0) - shift - 1; i >= 0; i--)
                        matrix[i + shift, j] = m[i, j];
                    for (int i = 0; i < shift; i++)
                        matrix[i, j] = 0;
                }
            }
            else Console.WriteLine("shift spot error");
            return matrix;
        }
        //----------------------------------------------------------------------------------------
        static double[,] shift_up(double[,] m, int shift)
        {
            double[,] matrix = new double[m.GetLength(0), m.GetLength(1)];
            if (shift < m.GetLength(0))
            {
                for (int j = 0; j < m.GetLength(1); j++)
                {
                    for (int i = 0; i < m.GetLength(0) - shift; i++)
                        matrix[i, j] = m[i + shift, j];
                    for (int i = m.GetLength(0) - shift; i < m.GetLength(0); i++)
                        matrix[i, j] = 0;
                }

            }
            else Console.WriteLine("shift spot error");
            return matrix;
        }
        //----------------------------------------------------------------------------------------
        #endregion
    }
}

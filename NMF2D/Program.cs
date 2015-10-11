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
        static double[,] XdiviXhat;
        static double[,] ones;


        static int tau; //時間ピッチ数
        static int fai; //周波数ピッチ数
        static int K;   //基底数
        static int I;
        static int J;

        static double min;

        static void Main(string[] args)
        {
            int divergence_flag = 0;           // 0:Euclid Norm  1:KL 
            int shift_flag = 0;                // 0:only f shift  1:f and time shift

            init(shift_flag);
            int itteration = 0;
            while (itteration < 9)
            {
                itteration++;
                Updates(divergence_flag);
                Console.WriteLine("itteration : " + itteration + " error = " + errorcalc(X, Xhat));

            }
            //CsvFileIO.CsvFileIO.WriteData("out_K.csv", T[0]);
            CsvFileIO.CsvFileIO.WriteData("reproducted.csv", reproduct(0));
        }
        //----------------------------------------------------------------------------
        static void init(int shift_flag)
        {
            string input_filename = @"C:\Users\優\Desktop\音素材\cq.csv";
            X = CsvFileIO.CsvFileIO.ReadData(input_filename);

            I = X.GetLength(0);
            J = X.GetLength(1);
            K = 3;
            tau = 7;
            fai = 12;

            XdiviXhat = new double[I, J];

            //T,Vの初期化
            switch (shift_flag)
            {
                // only f shift
                case 0:
                    tau = 1;
                    T.Add(initMatrix(I, K));
                    break;
                // f and time shift
                case 1:
                    for (int tt = 0; tt < tau; tt++)
                        T.Add(initMatrix(I, K));
                    break;

            }

            for (int ff = 0; ff < fai; ff++)
                V.Add(initMatrix(K, J));

            //onesの初期化
            ones = new double[I, J];
            for (int i = 0; i < I; i++)
                for (int j = 0; j < J; j++)
                    ones[i, j] = 1;


        }
        //----------------------------------------------------------------------------
        static void Updates(int divergence_flag)
        {
            updateXhat();
            updateT(divergence_flag);
            updateXhat();
            updateV(divergence_flag);
        }

        //----------------------------------------------------------------------------
        static void updateT(int divergence_flag)
        {
            double[,] A = new double[I, K];
            double[,] B = new double[I, K];

            //Euclid
            if (divergence_flag == 0)
            {
                for (int tt = 0; tt < tau; tt++)
                {
                    A.Clear();
                    B.Clear();
                    for (int ff = 0; ff < fai; ff++)
                    {
                        A = Mt.Add(A, Mt.Mul(shift_up(X, ff), shift_right(V[ff], tt).T()));
                        B = Mt.Add(B, Mt.Mul(shift_up(Xhat, ff), shift_right(V[ff], tt).T()));
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

            //KL
            else if (divergence_flag == 1)
            {
                for (int i = 0; i < I; i++)
                    for (int j = 0; j < J; j++)
                        if (Xhat[i, j] != 0)
                            XdiviXhat[i, j] = X[i, j] / Xhat[i, j];
                for (int tt = 0; tt < tau; tt++)
                {
                    A.Clear();
                    B.Clear();
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

        }
        //----------------------------------------------------------------------------
        static void updateV(int divergence_flag)
        {
            double[,] A = new double[K, J];
            double[,] B = new double[K, J];
            //Euclid
            if (divergence_flag == 0)
            {
                for (int ff = 0; ff < fai; ff++)
                {
                    A.Clear();
                    B.Clear();
                    for (int tt = 0; tt < tau; tt++)
                    {
                        A = Mt.Add(A, Mt.Mul(shift_down(T[tt], ff).T(), shift_left(X, tt)));
                        B = Mt.Add(B, Mt.Mul(shift_down(T[tt], ff).T(), shift_left(Xhat, tt)));
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

            //KL
            else if (divergence_flag == 1)
            {
                for (int i = 0; i < I; i++)
                    for (int j = 0; j < J; j++)
                        if (Xhat[i, j] != 0)
                            XdiviXhat[i, j] = X[i, j] / Xhat[i, j];

                for (int ff = 0; ff < fai; ff++)
                {
                    A.Clear();
                    B.Clear();
                    for (int tt = 0; tt < tau; tt++)
                    {
                        A = Mt.Add(A, Mt.Mul(shift_down(T[tt], ff).T(), shift_left(XdiviXhat, tt)));
                        B = Mt.Add(B, Mt.Mul(shift_down(T[tt], ff).T(), ones));
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
        static double[,] reproduct(int K)
        {
            double[,] A;
            double[,] rep = new double[I, J];
            for (int tt = 0; tt < tau; tt++)
            {
                for (int ff = 0; ff < fai; ff++)
                {
                    A = Mt.Mul(shift_down(getcolomn(T[tt], K), ff), shift_right(getrow(V[ff], K), tt));
                    rep = Mt.Add(rep, A);
                }
            }
            return rep;
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
        static double[,] getcolomn(double[,] m, int colomn_num)
        {
            double[,] output = new double[m.GetLength(0), 1];
            for (int i = 0; i < m.GetLength(0); i++)
                output[i, 0] = m[i, colomn_num];
            return output;
        }
        //----------------------------------------------------------------------------------------
        static double[,] getrow(double[,] m, int row_num)
        {
            double[,] output = new double[1, m.GetLength(1)];
            for (int i = 0; i < m.GetLength(1); i++)
                output[0, i] = m[row_num, i];
            return output;
        }
        //----------------------------------------------------------------------------------------
        static double errorcalc(double[,] truth, double[,] estimate)
        {
            int I = truth.GetLength(0);
            int J = truth.GetLength(1);
            double error = 0;
            for (int i = 0; i < I; i++)
            {
                for (int j = 0; j < J; j++)
                {
                    error += Math.Sqrt((truth[i, j] - estimate[i, j]) * (truth[i, j] - estimate[i, j]));
                }
            }
            return error / (I * J);
        }
        //----------------------------------------------------------------------------------------

        #endregion
    }
}

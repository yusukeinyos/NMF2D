using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using InoueLab;
using System.Numerics;

namespace MySignalProcessing
{
    public class SP
    {
        //-------------------------------------------------------------------------------------
        //実数値FFT（フル）
        //単独での使用
        public static Complex[] FullRealFastFourierTransform(double[] input)
        {
            input = zeroInserting(input); //0埋め
            Complex[] RFTout = Nm.RealFastFourierTransform(input);
            Complex[] output = new Complex[RFTout.Length * 2];
            for (int i = 0; i < RFTout.Length; i++)
            {
                output[i] = RFTout[i];
                output[RFTout.Length + i] = Complex.Conjugate(RFTout[RFTout.Length - 1 - i]);
            }
            return output;
        }
        //-------------------------------------------------------------------------------------
        //STFT（短時間フーリエ変換）
        public static Complex[,] STFT(double[] source_data, int window_size)
        {
            Complex[,] output;
            int source_length = source_data.Length; //入力データ長
            int T_step; //時間ステップ数
            int F_step; //周波数ステップ数
            int fft_start_point = 0; //fftをするスタート点
            double window_power = powerOfwindow(window_size); //窓関数のパワー WN

            T_step = (int)((source_length - window_size) / window_power) - 1;
            int n4 = 1 << (Mt.Log2Int(window_size) - 2);
            int n = n4 * 4, n2 = n4 * 2;
            F_step = n2 + 1;
            output = new Complex[F_step, T_step];

            for (int t = 0; t < T_step; t++)
            {
                double[] data = new double[window_size];
                fft_start_point = (int)window_power * t;
                for (int i = 0; i < window_size; i++)
                    data[i] = source_data[fft_start_point + i];
                data = windowing(data);

                Complex[] fft = Nm.RealFastFourierTransform(data);

                for (int f = 0; f < F_step; f++)
                    output[f, t] = fft[f];
            }

            return output;
        }
        //-------------------------------------------------------------------------------------
        //ISTFT(逆短時間フーリエ変換)
        public static double[] ISTFT(Complex[,] stft, int window_size)
        {
            int F_step = stft.GetLength(0);
            int T_step = stft.GetLength(1);
            double window_power = powerOfwindow(window_size); //窓関数のパワー WN
            double[] output = new double[(int)window_power * T_step + window_size]; //output
            double[] total_window = new double[output.Length]; //オーバーラップした窓関数の和

            double[] window = Hamming_Window(window_size);
            double[] ifft;
            Complex[] stft_t = new Complex[F_step];
            for (int i = 0; i < T_step; i++)
            {
                for (int f = 0; f < F_step; f++)
                    stft_t[f] = stft[f, i];
                Real_IFFT(stft_t, out ifft);

                for (int j = 0; j < window_size; j++)
                {
                    total_window[j + (int)window_power * i] += window[j];
                    output[j + (int)window_power * i] += ifft[j];
                }
            }
            for (int i = 0; i < output.Length; i++)
            {
                if (total_window[i] != 0)
                    output[i] /= total_window[i];
            }
            return output;
        }
        //-------------------------------------------------------------------------------------
        //実数値IFFT（高速フーリエ逆変換）
        //複素共役にする処理もNm.RealFastFourierTransform(Complex[])に含まれてる!!
        //入力はナイキスト周波数に対応するスペクトル（半分だけ）でいい!!
        public static void Real_IFFT(double[] Re, double[] Im, out double[] output)
        {
            int size = Re.Length;
            Complex[] comp = new Complex[size];

            for (int i = 0; i < size; i++)
                comp[i] = new Complex(Re[i], Im[i]);
            output = Nm.RealFastFourierTransform(comp);

        }
        static void Real_IFFT(Complex[] comp, out double[] output)
        {
            output = Nm.RealFastFourierTransform(comp);
        }
        //-------------------------------------------------------------------------------------
        //スペクトログラム　⇒　対数スペクトログラム 10*log_10(STFT^2) スケール変換あり
        public static double[,] logSpectrogram(double[,] spectrogram, out double min)
        {
            min = 1;
            for (int i = 0; i < spectrogram.GetLength(0); i++)
                for (int j = 0; j < spectrogram.GetLength(1); j++)
                {
                    if (min > spectrogram[i, j] && spectrogram[i, j] != 0)
                        min = spectrogram[i, j];
                }
            for (int i = 0; i < spectrogram.GetLength(0); i++)
                for (int j = 0; j < spectrogram.GetLength(1); j++)
                {
                    if (spectrogram[i, j] != 0)
                        spectrogram[i, j] = 1 * Math.Log10(spectrogram[i, j]);
                    else
                        spectrogram[i, j] = 1 * Math.Log10(min);
                    spectrogram[i, j] -= 1 * Math.Log10(min);
                }
            return spectrogram;
        }
        //---------------------------------------------------------------------------------
        //対数スペクトログラム ⇒　スペクトログラム
        public static double[,] Inv_logSpectrogram(double[,] spectrogram, double min)
        {
            for (int i = 0; i < spectrogram.GetLength(0); i++)
                for (int j = 0; j < spectrogram.GetLength(1); j++)
                {
                    spectrogram[i, j] += Math.Log10(min);
                    spectrogram[i, j] = Math.Pow(10, spectrogram[i, j] / 1);
                }
            return spectrogram;
        }
        //---------------------------------------------------------------------------------
        #region Parts
        //-------------------------------------------------------------
        //0埋め
        private static double[] zeroInserting(double[] data)
        {
            int size = data.Length;
            int fixed_size = 1 << Mt.Log2Int(size) + 1;
            if (fixed_size > size)
            {
                double[] outdata = new double[fixed_size];
                for (int i = 0; i < size; i++)
                    outdata[i] = data[i];

                return outdata;
            }
            else
                return data;
        }
        //-------------------------------------------------------------
        //窓関数を掛ける
        static double[] windowing(double[] input)
        {
            int window_size = input.Length; //窓幅
            double[] window = Hamming_Window(window_size); //窓関数

            for (int i = 0; i < window_size; i++)
                input[i] *= window[i];

            return input;
        }
        //-------------------------------------------------------------
        //窓関数
        static double[] Hamming_Window(int window_size)
        {
            double[] output = new double[window_size];
            for (int i = 0; i < window_size; i++)
            {
                double h = 2 * Math.PI / (window_size - 1);
                output[i] = 0.54 - 0.46 * Math.Cos(h * i);
            }
            return output;
        }
        //-------------------------------------------------------------
        //WNを返す
        static double powerOfwindow(int window_size)
        {
            double[] window = Hamming_Window(window_size);
            double sum = 0;
            for (int i = 0; i < window_size; i++)
                sum += window[i] * window[i];
            return sum;
        }
        //-------------------------------------------------------------
        //複素行列　⇒　パワー行列
        public static double[,] CmatToPowermat(Complex[,] c)
        {
            double[,] d = new double[c.GetLength(0), c.GetLength(1)];
            for (int i = 0; i < c.GetLength(0); i++)
                for (int j = 0; j < c.GetLength(1); j++)
                    d[i, j] = c[i, j].Magnitude * c[i, j].Magnitude;
            return d;
        }
        //-------------------------------------------------------------
        //パワー行列d＋複素行列cの位相成分　⇒　複素行列
        public static Complex[,] PowermatToCmat(double[,] d, Complex[,] c)
        {
            Complex[,] output = new Complex[c.GetLength(0), c.GetLength(1)];
            for (int i = 0; i < c.GetLength(0); i++)
                for (int j = 0; j < c.GetLength(1); j++)
                    output[i, j] = new Complex(Math.Sqrt(d[i, j]) * Math.Cos(c[i, j].Phase), Math.Sqrt(d[i, j]) * Math.Sin(c[i, j].Phase));
            return output;
        }
        //-------------------------------------------------------------
        //モノラルデータ⇒ステレオデータ（ジャグ配列）
        public static void arrayTojag(double[] d, out double[][] output)
        {
            output = new double[2][];
            output[0] = (double[])d.Clone();
            output[1] = (double[])d.Clone();
        }
        //-------------------------------------------------------------
        # endregion
    }
}

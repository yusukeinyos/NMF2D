using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Numerics;  // projectの参照設定でSystem.Numericsを追加
using System.Text;
using InoueLab;

namespace System.Linq
{
    public static partial class EnumerableEx
    {
        public static BigInteger Sum<TSource>(this IEnumerable<TSource> source, Func<TSource, BigInteger> selector) { return source.Select(element => selector(element)).Sum(); }
        public static BigInteger Sum(this IEnumerable<BigInteger> source)
        {
            BigInteger sum = 0;
            foreach (var element in source) checked { sum += element; }
            return sum;
        }
        public static Complex Sum<TSource>(this IEnumerable<TSource> source, Func<TSource, Complex> selector) { return source.Select(element => selector(element)).Sum(); }
        public static Complex Sum(this IEnumerable<Complex> source)
        {
            Complex sum = 0;
            foreach (var element in source) checked { sum += element; }
            return sum;
        }
        public static Complex Average<TSource>(this IEnumerable<TSource> source, Func<TSource, Complex> selector) { return source.Select(element => selector(element)).Average(); }
        public static Complex Average(this IEnumerable<Complex> source)
        {
            Complex sum = 0;
            int count = 0;
            foreach (var element in source)
            {
                checked { sum += element; }
                count++;
            }
            return sum / count;
        }
    }
}

namespace InoueLab
{
    #region Spline classes
    public class Spline3
    {
        double[] DataX, DataY, DataZ;
        public Spline3() { }
        public Spline3(IEnumerable<double> dataX, IEnumerable<double> dataY) { Set(dataX, dataY); }

        int FindSection(double value)
        {
            int i = 0;
            for (int j = DataX.Length - 1; i < j; )
            {
                int k = (i + j) / 2;
                if (DataX[k] < value) i = k + 1; else j = k;
            }
            if (i > 0) --i;
            return i;
        }
        public double Interpolate(double value)
        {
            int i = FindSection(value);
            double h = DataX[i + 1] - DataX[i];
            double d = value - DataX[i];
            return (((DataZ[i + 1] - DataZ[i]) * d / h + DataZ[i] * 3) * d + ((DataY[i + 1] - DataY[i]) / h - (DataZ[i] * 2 + DataZ[i + 1]) * h)) * d + DataY[i];
        }
        public double CalcGrad(double value)
        {
            int i = FindSection(value);
            double h = DataX[i + 1] - DataX[i];
            double d = value - DataX[i];
            return ((DataZ[i + 1] - DataZ[i]) * d / h * 3 + DataZ[i] * 6) * d + ((DataY[i + 1] - DataY[i]) / h - (DataZ[i] * 2 + DataZ[i + 1]) * h);
        }

        public void Set(IEnumerable<double> dataX, IEnumerable<double> dataY)
        {
            if (dataY == null) ThrowException.ArgumentException("dataY");
            DataY = dataY.ToArray();
            int N = DataY.Length;
            if (N == 0) ThrowException.ArgumentException("dataY");
            if (N == 1) DataY = new double[] { DataY[0], DataY[0], DataY[0] };
            if (N == 2) DataY = new double[] { DataY[0], (DataY[0] + DataY[1]) / 2, DataY[1] };

            if (dataX != null)
            {
                DataX = dataX.ToArray();
                if (DataX.Length == 1) DataX = new double[] { DataX[0], DataX[0], DataX[0] };
                if (DataX.Length == 2) DataX = new double[] { DataX[0], (DataX[0] + DataX[1]) / 2, DataX[1] };
            }
            if (dataX == null || DataX.Length == 0)
            {
                DataX = new double[N];
                DataX[0] = 0;
                for (int i = 1; i < N; i++) DataX[i] = DataX[i - 1] + Math.Max(Math.Abs(DataY[i] - DataY[i - 1]), 1e-10);
                for (int i = 1; i < N; i++) DataX[i] /= DataX[N - 1];
            }
            if (DataX.Length != N) ThrowException.ArgumentOutOfRangeException("dataX, dataY");

            DataZ = new double[N];
            double[] h = new double[N];
            double[] d = new double[N];
            DataZ[0] = DataZ[N - 1] = 0;
            for (int i = 0; i < N - 1; i++)
            {
                h[i] = DataX[i + 1] - DataX[i];
                d[i + 1] = (DataY[i + 1] - DataY[i]) / h[i];
            }
            DataZ[1] = d[2] - d[1] - h[0] * DataZ[0];
            d[1] = 2 * (DataX[2] - DataX[0]);
            for (int i = 1; i < N - 2; i++)
            {
                double t = h[i] / d[i];
                DataZ[i + 1] = d[i + 2] - d[i + 1] - DataZ[i] * t;
                d[i + 1] = 2 * (DataX[i + 2] - DataX[i]) - h[i] * t;
            }
            DataZ[N - 2] -= h[N - 2] * DataZ[N - 1];
            for (int i = N - 2; i > 0; i--)
                DataZ[i] = (DataZ[i] - h[i] * DataZ[i + 1]) / d[i];
        }
    }

    public class Spline1
    {
        double[] DataX, DataY;
        public Spline1() { }
        public Spline1(IEnumerable<double> dataX, IEnumerable<double> dataY) { Set(dataX, dataY); }

        int FindSection(double value)
        {
            int i = 0;
            for (int j = DataX.Length - 1; i < j; )
            {
                int k = (i + j) / 2;
                if (DataX[k] < value) i = k + 1; else j = k;
            }
            if (i > 0) --i;
            return i;
        }
        public double Interpolate(double value)
        {
            int i = FindSection(value);
            return (DataY[i + 1] - DataY[i]) / (DataX[i + 1] - DataX[i]) * (value - DataX[i]) + DataY[i];
        }
        public void Set(IEnumerable<double> dataX, IEnumerable<double> dataY)
        {
            if (dataX == null || dataY == null) ThrowException.ArgumentException("dataX, dataY");
            DataX = dataX.ToArray();
            DataY = dataY.ToArray();
            if (DataX.Length == 0 || DataY.Length == 0 || DataX.Length != DataY.Length) ThrowException.ArgumentException("dataX, dataY");
        }
    }
    #endregion

    #region Numerical Integrator
    public class Integrator
    {
        const int LoopMax = 14, K = 7;
        const double EPS = 1e-8;

        public bool Continue = true;
        public double left, right, width;
        public double Value { get { return _Value; } }
        public double Answer { get { return _Answer; } }
        public double Error { get { return _Error; } }

        double _Value;
        double _Answer;
        double _Error;
        int counter;
        int binN;
        double sum = 0;
        List<double> dataX = new List<double>();
        List<double> dataY = new List<double>();

        public Integrator() { }
        public Integrator(double left, double right) { Init(left, right); }

        public void Init(double left, double right)
        {
            this.left = left; this.right = right; width = right - left;
            _Value = left + width / 2;
            counter = 1;
        }

        static Double2 polint(List<double> dataX, List<double> dataY, int start, int count)
        {
            double[] y0 = Enumerable.Range(start, count).Select(i => dataY[start + i]).ToArray();
            double[] y1 = y0.CloneX();
            int j = count;
            double y = y1[--j], dy = 0;
            for (int m = 1; m < count; m++)
            {
                for (int i = 0; i < count - m; i++)
                {
                    double x0 = dataX[start + i];
                    double x1 = dataX[start + i + m];
                    double k = (y0[i + 1] - y1[i]) / (x0 - x1);
                    y0[i] = x0 * k;
                    y1[i] = x1 * k;
                }
                dy = (2 * j < count - m) ? y0[j] : y1[--j];
                y += dy;
            }
            return new Double2(y, dy);
        }
        void Input_(double value)
        {
            dataX.Add(Math.Pow(9, -dataX.Count));
            dataY.Add(dataY.Count == 0 ? width * sum : (dataY[dataY.Count - 1] + width * sum / (binN / 6)) / 3);

            if (dataY.Count >= K)
            {
                Double2 r = polint(dataX, dataY, dataY.Count - K, K);
                _Answer = r.X; _Error = r.Y;
                if (Continue)
                {
                    if (Math.Abs(_Error) <= EPS * Math.Abs(_Answer) ||
                        _Error * _Error <= EPS * EPS * Math.Abs(_Answer) ||
                        dataY.Count == LoopMax)
                        Continue = false;
                }
            }

            binN = (int)Math.Pow(3, dataY.Count) * 2;
            counter = binN / 3;
            sum = 0;
        }
        public bool Input(double value)
        {
            sum += value;
            if (--counter == 0) Input_(value);
            _Value = left + width * (counter / 2 * 6 + ((counter & 1) * 4 + 1)) / binN;
            return Continue;
        }
        public void Break() { Continue = false; }
    }
    #endregion

    public static partial class Ex
    {
        public static BigInteger Product<TSource>(this IEnumerable<TSource> source, Func<TSource, BigInteger> selector) { return source.Select(element => selector(element)).Product(); }
        public static BigInteger Product(this IEnumerable<BigInteger> source)
        {
            BigInteger product = 1;
            foreach (var element in source) checked { product *= element; }
            return product;
        }
        public static Complex Product<TSource>(this IEnumerable<TSource> source, Func<TSource, Complex> selector) { return source.Select(element => selector(element)).Product(); }
        public static Complex Product(this IEnumerable<Complex> source)
        {
            Complex product = 1;
            foreach (var element in source) checked { product *= element; }
            return product;
        }

        public static void Sort<T>(this List<T> list, Func<T, BigInteger> selector) { list.Sort((x, y) => selector(x).CompareTo(selector(y))); }
        public static void Sort<T>(this List<T> list, Func<T, BigInteger> selector, bool ad) { if (ad) list.Sort((x, y) => selector(x).CompareTo(selector(y))); else list.Sort((x, y) => selector(y).CompareTo(selector(x))); }
        public static int[] SortIndex<T>(this IList<T> list, Func<T, BigInteger> selector) { return SortIndex(list, (x, y) => selector(x).CompareTo(selector(y))); }
        public static int[] SortIndex<T>(this IList<T> list, Func<T, BigInteger> selector, bool ad) { return ad ? SortIndex(list, (x, y) => selector(x).CompareTo(selector(y))) : SortIndex(list, (x, y) => selector(y).CompareTo(selector(x))); }
    }

    public static partial class Mt
    {
        public static BigInteger MultinomialInteger(IEnumerable<int> source)
        {
            int total = 0;
            BigInteger product = 1;
            foreach (int element in source)
            {
                if (element < 0) ThrowException.ArgumentException("element");
                total += element;
                product *= FactorialInteger(element);
            }
            return FactorialInteger(total) / product;
        }

        static List<BigInteger> FactorialIntegerBuffer = new List<BigInteger>() { 1 };
        static BigInteger FactorialInteger_(int value)
        {
            BigInteger product = FactorialIntegerBuffer.Last();
            for (int i = FactorialIntegerBuffer.Count; i <= value; i++)
            {
                product *= i;
                FactorialIntegerBuffer.Add(product);
            }
            return product;
        }
        public static BigInteger FactorialInteger(int value)
        {
            if (value < 0) ThrowException.ArgumentOutOfRangeException("value");
            return value < FactorialIntegerBuffer.Count ? FactorialIntegerBuffer[value] : FactorialInteger_(value);
        }

        public static BigInteger GreatestCommonDivisor(BigInteger val1, BigInteger val2)
        {
            while (true)
            {
                if (val1 < val2) Mt.Swap(ref val1, ref val2);
                var z = BigInteger.Remainder(val1, val2);
                if (z == 0) break;
                val1 = z;
            }
            return val2;
        }
    }

    public static partial class Nm
    {
        #region Signal processing functions
        public enum DataWindowType { Box, Hanning, Hamming, Blackman, Parzen, Welch, NormalDistribution = -1 }
        public static double[] GetDataWindow(DataWindowType type, int size)
        {
            double[] Table = new double[size];
            Table = new double[size];

            double h = 2 * Math.PI / (size - 1);
            switch (type)
            {
                case DataWindowType.Box:
                    for (int i = size; --i >= 0; ) Table[i] = 1;
                    break;
                case DataWindowType.Hanning:
                    for (int i = size; --i >= 0; ) Table[i] = 0.50 - 0.50 * Math.Cos(h * i);
                    break;
                case DataWindowType.Hamming:
                    for (int i = size; --i >= 0; ) Table[i] = 0.54 - 0.46 * Math.Cos(h * i);
                    break;
                case DataWindowType.Blackman:
                    for (int i = size; --i >= 0; ) Table[i] = 0.42 - 0.50 * Math.Cos(h * i) + 0.08 * Math.Cos(2 * h * i);
                    break;
                case DataWindowType.Parzen:
                    for (int i = size; --i >= 0; ) Table[i] = 1.0 - Math.Abs((i * 2 - (size - 1)) / (double)(size + 1));
                    break;
                case DataWindowType.Welch:
                    for (int i = size; --i >= 0; ) Table[i] = 1.0 - Mt.Sq((i * 2 - (size - 1)) / (double)(size + 1));
                    break;
            }
            //double c = Math.Sqrt(n / Table.Sum(x => Sq(x)));
            //for (int i = n; --i >= 0; ) Table[i] *= c;
            return Table;
        }

        public static void Windowing(double[] data, DataWindowType type)
        {
            double[] window = GetDataWindow(type, data.Length);
            for (int i = data.Length; --i >= 0; ) data[i] *= window[i];
        }


        public class TriangularTable
        {
            double[] Table;
            int size, size4, size2, size3;

            TriangularTable(int size)
            {
                if (size <= 0 || size % 4 != 0) ThrowException.ArgumentOutOfRangeException("size");
                this.size4 = size / 4;
                this.size2 = size4 * 2;
                this.size3 = size4 * 3;
                this.size = size4 * 4;
                Table = new double[size4 + 1];

                Table[size4] = 1;
                double t = 2 * Math.PI / size;
                for (int i = size4; --i > 0; )
                    Table[i] = Math.Sin(t * i);
                Table[0] = 0;
            }

            public double Cos(int index)
            {
                if (index < 0) index = index % size + size; else if (index >= size) index %= size;
                if (index < size4) return Table[size4 - index];
                if (index < size2) return -Table[index - size4];
                if (index < size3) return -Table[size3 - index];
                return Table[index - size3];
                //return (i < n4 || i >= n3 ? 1 : -1) * Table[i < n2 ? (i < n4 ? n4 - i : i - n4) : (i < n3 ? n3 - i : i - n3)];
            }
            public double Sin(int index)
            {
                if (index < 0) index = index % size + size; else if (index >= size) index %= size;
                if (index < size4) return Table[index];
                if (index < size2) return Table[size2 - index];
                if (index < size3) return -Table[index - size2];
                return -Table[size - index];
                //return (i < n2 ? 1 : -1) * Table[i <= n2 ? (i < n4 ? i : n2 - i) : (i < n3 ? i - n2 : n - i)];
            }
            public Complex Complex(int index) { return new Complex(Cos(index), Sin(index)); }

            static List<TriangularTable> Stocks = new List<TriangularTable>();
            public static TriangularTable Get(int size)
            {
                for (int i = Stocks.Count; --i >= 0; )
                {
                    var s = Stocks[i];
                    if (s.size == size)
                    {
                        if (i < Stocks.Count - 1)
                        {
                            Stocks.RemoveAt(i);
                            Stocks.Add(s);
                        }
                        return s;
                    }
                }
                {
                    var s = new TriangularTable(size);
                    if (Stocks.Count >= 8) Stocks.RemoveAt(0);
                    Stocks.Add(s);
                    return s;
                }
            }
        }

        static Complex[] FastFourierTransform_(Complex[] src, bool sw)
        {
            int l = Mt.Log2Int(src.Length);
            int n = 1 << l;
            Complex[] Data = new Complex[n];
            for (int j = -(n >> 1), i = 0; i < n; i++)
            {
                int k = n;
                while ((k >>= 1) <= j) j -= k;
                j += k;
                Data[i] = src[j];
            }

            int sign = sw ? -1 : 1;
            TriangularTable Triangle = TriangularTable.Get(n);
            for (int m = 0; m < l; m++)
            {
                int step = 1 << m;
                int rotw = sign * (1 << (l - 1 - m));
                for (int k = 0; k < step; k++)
                {
                    Complex u = Triangle.Complex(rotw * k);  // (2 * Math.PI / n) * rotw * k
                    for (int i = k; i < n; i += step * 2)
                    {
                        Complex t = Data[i + step] * u;
                        Data[i + step] = Data[i] - t;
                        Data[i] += t;
                    }
                }
            }
            return Data;
        }

        public static Complex[] FastFourierTransform(Complex[] src, bool sw, double amplitude = 1.0)
        {
            Complex[] Data = FastFourierTransform_(src, sw);
            if (!sw) amplitude /= Data.Length;
            LetMul(Data, amplitude);
            return Data;
        }

        public static Complex[] RealFastFourierTransform(double[] src, double amplitude = 1.0)
        {
            int n4 = 1 << (Mt.Log2Int(src.Length) - 2);
            int n = n4 * 4, n2 = n4 * 2, n3 = n4 * 3;
            int n2mask = n2 - 1;

            Complex[] Data = new Complex[n2 + 1];
            for (int i = n2; --i >= 0; ) Data[i] = new Complex(src[i * 2], src[i * 2 + 1]);
            Complex[] Dst = FastFourierTransform_(Data, true);
            TriangularTable Triangle = TriangularTable.Get(n);
            for (int i = n4; i >= 0; --i)
            {
                int j = n2 - i;
                Complex g1 = Dst[i];
                Complex g2 = Complex.Conjugate(Dst[j & n2mask]);
                Complex h1 = (g1 + g2);
                Complex h2 = (g1 - g2) * Triangle.Complex(n3 - i);
                Data[i] = (h1 + h2);
                Data[j] = Complex.Conjugate(h1 - h2);
            }
            LetMul(Data, amplitude / 2);
            return Data;
        }
        public static double[] RealFastFourierTransform(Complex[] src, double amplitude = 1.0)
        {
            int n4 = 1 << (Mt.Log2Int(src.Length) - 1);
            int n = n4 * 4, n2 = n4 * 2, n3 = n4 * 3;
            int n2mask = n2 - 1;

            Complex[] Data = new Complex[n2];
            TriangularTable Triangle = TriangularTable.Get(n);
            for (int i = n4; i >= 0; --i)
            {
                int j = n2 - i;
                Complex g1 = src[i];
                Complex g2 = Complex.Conjugate(src[j]);
                Complex h1 = (g1 + g2);
                Complex h2 = (g1 - g2) * Triangle.Complex(n4 + i);
                Data[i] = (h1 + h2);
                Data[j & n2mask] = Complex.Conjugate(h1 - h2);
            }
            Data = FastFourierTransform_(Data, false);
            LetMul(Data, amplitude / n);
            double[] Dst = New.Array(n, i => (i & 1) == 0 ? Data[i / 2].Real : Data[i / 2].Imaginary);
            return Dst;
        }
        static void LetMul(Complex[] left, double right)
        {
            if (right == 1) return;
            for (int i = left.Length; --i >= 0; )
                left[i] = new Complex(left[i].Real * right, left[i].Imaginary * right);
        }

        // Power Spectral Densityを求めたい。nの増加に伴いそれぞれの周波数カラムの幅(Δf)は狭くなるため、結果のData値は増加する
        // ΣResult   != 単位時間辺りのpower
        // ∫Result df = 単位時間辺りのpower
        // y=a*sin(x) の結果は、a*a*n/2 振幅はa/2でパワーはa*a/4、one-sidedなので2倍してa*a/2、nは上記
        // これは (1/n)*Σ(y*y) の結果と一致する
        // y=a の結果は a*a*n
        public static double[] PowerSpectrumFFT(double[] data, double amplitude)
        {
            Complex[] Freq = RealFastFourierTransform(data, Math.Sqrt(2 * amplitude / data.Length));
            double[] Powr = new double[Freq.Length];
            for (int i = Powr.Length; --i >= 0; ) Powr[i] = Mt.Sq(Freq[i].Real) + Mt.Sq(Freq[i].Imaginary);
            Powr[0] /= 2;
            Powr[Powr.Length - 1] /= 2;
            return Powr;
        }

        public static Complex[] PowerPhaseSpectrumFFT(double[] data, double amplitude)
        {
            Complex[] Freq = RealFastFourierTransform(data, Math.Sqrt(2 * amplitude / data.Length));
            Complex[] Powr = new Complex[Freq.Length];
            for (int i = Powr.Length; --i >= 0; ) Powr[i] = new Complex(Mt.Sq(Freq[i].Real) + Mt.Sq(Freq[i].Imaginary), Freq[i].Phase);
            var a = Powr[0]; Powr[0] = new Complex(a.Real / 2, a.Imaginary);
            var b = Powr[Powr.Length - 1]; Powr[Powr.Length - 1] = new Complex(b.Real / 2, b.Imaginary);
            return Powr;
        }
        #endregion

        //---------------------------------------------------------------------------
        #region Numerical optimization functions
        public static double[] ArgminSimplex(double[] arguments, double[] delta, Func<double[], double> func, double tolerance)
        {
            int EstimationCountMax = 5000;
            int EstimationCount = 0;
            int D = arguments.Length, D1 = D + 1;

            double[][] Args = New.Array(D + 1, j => New.Array(D, i => arguments[i] + (i == j ? delta[i] : 0.0)));
            double[] ArgSum = New.Array(D, i => Ex.Range(D + 1).Sum(j => Args[j][i]));
            double[] ArgTry = new double[D];
            double[] Values = Args.Select(arg => func(arg)).ToArray();

            Func<int, double, double> tryfunc = (int i, double fac) =>
            {
                double[] arg = Args[i];
                double fac1 = (1 - fac) / D;
                double fac2 = fac1 - fac;
                for (int j = D; --j >= 0; )
                    ArgTry[j] = ArgSum[j] * fac1 - arg[j] * fac2;

                double y = func(ArgTry);
                if (Values[i] > y)
                {
                    Values[i] = y;
                    for (int j = D; --j >= 0; )
                    {
                        ArgSum[j] += ArgTry[j] - arg[j];
                        arg[j] = ArgTry[j];
                    }
                }
                return y;
            };

            while (true)
            {
                int iMin = 0;
                int iMax = Values[0] > Values[1] ? 0 : 1;
                int iMax2 = 1 - iMax;
                for (int i = 0; i < D1; i++)
                {
                    double v = Values[i];
                    if (v <= Values[iMin]) iMin = i;
                    if (v > Values[iMax]) { iMax2 = iMax; iMax = i; }
                    else if (v > Values[iMax2] && i != iMax) iMax2 = i;
                }

                double d = 2.0 * Math.Abs(Values[iMax] - Values[iMin])
                    / (Math.Abs(Values[iMax]) + Math.Abs(Values[iMin]) + Mt.DoubleEps);
                if (d < tolerance) break;
                if (EstimationCount >= EstimationCountMax)
                {
                    ThrowException.WriteLine("ArgminSimplex: too many function estimations.");
                    break;
                }

                EstimationCount += 2;
                double vTry = tryfunc(iMax, -1);
                if (vTry <= Values[iMin]) { tryfunc(iMax, 2); continue; }
                if (vTry < Values[iMax2]) { EstimationCount--; continue; }

                double vTryPrev = Values[iMax];
                vTry = tryfunc(iMax, 0.5);
                if (vTry < vTryPrev) continue;

                EstimationCount += D;
                for (int i = 0; i < D1; i++)
                {
                    if (i == iMin) continue;
                    for (int j = 0; j < D; j++)
                        Args[i][j] = (Args[i][j] + Args[iMin][j]) / 2;
                    Values[i] = func(Args[i]);
                }
                for (int i = D; --i >= 0; )
                    ArgSum[i] = Ex.Range(D + 1).Sum(j => Args[j][i]);
            }
            return New.Array(D, i => Ex.Range(D1).Average(j => Args[j][i]));
        }
        #endregion

    }

}
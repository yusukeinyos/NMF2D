using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace InoueLab
{
    #region Error class
    public static class ThrowException
    {
        public static void WriteLine(string message)
        {
            Console.Error.WriteLine(message);
        }
        public static void ArgumentException(string message)
        {
            throw new ArgumentException(message);
        }
        public static void ArgumentNullException(string paramName)
        {
            throw new ArgumentNullException(paramName);
        }
        public static void ArgumentOutOfRangeException(string paramName)
        {
            throw new ArgumentOutOfRangeException(paramName);
        }
        public static void InvalidOperationException(string message)
        {
            throw new InvalidOperationException(message);
        }
        public static void NoElements()
        {
            throw new InvalidOperationException("no elements");
        }
        public static void SizeMismatch()
        {
            throw new InvalidOperationException("size mismatch");
        }

    }
    #endregion

    #region Concatenated classes
    public struct V2<T0, T1>
    {
        public T0 v0;
        public T1 v1;
        public V2(T0 v0, T1 v1) { this.v0 = v0; this.v1 = v1; }
        public override string ToString() { return string.Format("{0}, {1}", v0, v1); }
        //public override bool Equals(object obj)
        //{
        //    if (obj is V2<T0, T1>) {
        //        V2<T0, T1> right = (V2<T0, T1>)obj;
        //        if (v0.Equals(v0) && v1.Equals(right.v1)) return true;
        //    }
        //    return false;
        //}
        //public static bool operator ==(V2<T0, T1> left, V2<T0, T1> right) { return left.v0.Equals(right.v0) && left.v1.Equals(right.v1); }
    }

    public struct V3<T0, T1, T2>
    {
        public T0 v0;
        public T1 v1;
        public T2 v2;
        public V3(T0 v0, T1 v1, T2 v2) { this.v0 = v0; this.v1 = v1; this.v2 = v2; }
        public override string ToString() { return string.Format("{0}, {1}, {2}", v0, v1, v2); }
    }

    public struct V4<T0, T1, T2, T3>
    {
        public T0 v0;
        public T1 v1;
        public T2 v2;
        public T3 v3;
        public V4(T0 v0, T1 v1, T2 v2, T3 v3) { this.v0 = v0; this.v1 = v1; this.v2 = v2; this.v3 = v3; }
        public override string ToString() { return string.Format("{0}, {1}, {2}, {3}", v0, v1, v2, v3); }
    }

    public struct Array1<T> : IList<T>
    {
        public T v0;
        public Array1(T v0)
        {
            this.v0 = v0;
        }

        public int IndexOf(T item)
        {
            if (item.Equals(v0)) return 0;
            return -1;
        }
        public void Insert(int index, T item) { throw new NotImplementedException(); }
        public void RemoveAt(int index) { throw new NotImplementedException(); }
        public T this[int index]
        {
            get
            {
                switch (index)
                {
                    case 0: return v0;
                }
                ThrowException.ArgumentException("index"); return default(T);
            }
            set
            {
                switch (index)
                {
                    case 0: v0 = value; return;
                }
                ThrowException.ArgumentException("index");
            }
        }
        public void Add(T item) { throw new NotImplementedException(); }
        public void Clear() { throw new NotImplementedException(); }
        public bool Contains(T item) { return IndexOf(item) != -1; }
        public void CopyTo(T[] array, int arrayIndex)
        {
            array[arrayIndex] = v0;
        }
        public int Count { get { return 1; } }
        public bool IsReadOnly { get { return false; } }
        public bool Remove(T item) { throw new NotImplementedException(); }
        public IEnumerator<T> GetEnumerator()
        {
            yield return v0;
        }
        IEnumerator IEnumerable.GetEnumerator() { return GetEnumerator(); }
    }

    public struct Array2<T> : IList<T>
    {
        public T v0;
        public T v1;
        public Array2(T v0, T v1)
        {
            this.v0 = v0;
            this.v1 = v1;
        }

        public int IndexOf(T item)
        {
            if (item.Equals(v0)) return 0;
            if (item.Equals(v1)) return 1;
            return -1;
        }
        public void Insert(int index, T item) { throw new NotImplementedException(); }
        public void RemoveAt(int index) { throw new NotImplementedException(); }
        public T this[int index]
        {
            get
            {
                switch (index)
                {
                    case 0: return v0;
                    case 1: return v1;
                }
                ThrowException.ArgumentException("index"); return default(T);
            }
            set
            {
                switch (index)
                {
                    case 0: v0 = value; return;
                    case 1: v1 = value; return;
                }
                ThrowException.ArgumentException("index");
            }
        }
        public void Add(T item) { throw new NotImplementedException(); }
        public void Clear() { throw new NotImplementedException(); }
        public bool Contains(T item) { return IndexOf(item) != -1; }
        public void CopyTo(T[] array, int arrayIndex)
        {
            array[arrayIndex] = v0;
            array[arrayIndex + 1] = v1;
        }
        public int Count { get { return 2; } }
        public bool IsReadOnly { get { return false; } }
        public bool Remove(T item) { throw new NotImplementedException(); }
        public IEnumerator<T> GetEnumerator()
        {
            yield return v0;
            yield return v1;
        }
        IEnumerator IEnumerable.GetEnumerator() { return GetEnumerator(); }
    }

    public struct Array3<T> : IList<T>
    {
        public T v0;
        public T v1;
        public T v2;
        public Array3(T v0, T v1, T v2)
        {
            this.v0 = v0;
            this.v1 = v1;
            this.v2 = v2;
        }

        public int IndexOf(T item)
        {
            if (item.Equals(v0)) return 0;
            if (item.Equals(v1)) return 1;
            if (item.Equals(v2)) return 2;
            return -1;
        }
        public void Insert(int index, T item) { throw new NotImplementedException(); }
        public void RemoveAt(int index) { throw new NotImplementedException(); }
        public T this[int index]
        {
            get
            {
                switch (index)
                {
                    case 0: return v0;
                    case 1: return v1;
                    case 2: return v2;
                }
                ThrowException.ArgumentException("index"); return default(T);
            }
            set
            {
                switch (index)
                {
                    case 0: v0 = value; return;
                    case 1: v1 = value; return;
                    case 2: v2 = value; return;
                }
                ThrowException.ArgumentException("index");
            }
        }
        public void Add(T item) { throw new NotImplementedException(); }
        public void Clear() { throw new NotImplementedException(); }
        public bool Contains(T item) { return IndexOf(item) != -1; }
        public void CopyTo(T[] array, int arrayIndex)
        {
            array[arrayIndex] = v0;
            array[arrayIndex + 1] = v1;
            array[arrayIndex + 2] = v2;
        }
        public int Count { get { return 2; } }
        public bool IsReadOnly { get { return false; } }
        public bool Remove(T item) { throw new NotImplementedException(); }
        public IEnumerator<T> GetEnumerator()
        {
            yield return v0;
            yield return v1;
            yield return v2;
        }
        IEnumerator IEnumerable.GetEnumerator() { return GetEnumerator(); }
    }

    public struct Array4<T> : IList<T>
    {
        public T v0;
        public T v1;
        public T v2;
        public T v3;
        public Array4(T v0, T v1, T v2, T v3)
        {
            this.v0 = v0;
            this.v1 = v1;
            this.v2 = v2;
            this.v3 = v3;
        }

        public int IndexOf(T item)
        {
            if (item.Equals(v0)) return 0;
            if (item.Equals(v1)) return 1;
            if (item.Equals(v2)) return 2;
            if (item.Equals(v3)) return 3;
            return -1;
        }
        public void Insert(int index, T item) { throw new NotImplementedException(); }
        public void RemoveAt(int index) { throw new NotImplementedException(); }
        public T this[int index]
        {
            get
            {
                switch (index)
                {
                    case 0: return v0;
                    case 1: return v1;
                    case 2: return v2;
                    case 3: return v3;
                }
                ThrowException.ArgumentException("index"); return default(T);
            }
            set
            {
                switch (index)
                {
                    case 0: v0 = value; return;
                    case 1: v1 = value; return;
                    case 2: v2 = value; return;
                    case 3: v3 = value; return;
                }
                ThrowException.ArgumentException("index");
            }
        }
        public void Add(T item) { throw new NotImplementedException(); }
        public void Clear() { throw new NotImplementedException(); }
        public bool Contains(T item) { return IndexOf(item) != -1; }
        public void CopyTo(T[] array, int arrayIndex)
        {
            array[arrayIndex] = v0;
            array[arrayIndex + 1] = v1;
            array[arrayIndex + 2] = v2;
            array[arrayIndex + 3] = v3;
        }
        public int Count { get { return 2; } }
        public bool IsReadOnly { get { return false; } }
        public bool Remove(T item) { throw new NotImplementedException(); }
        public IEnumerator<T> GetEnumerator()
        {
            yield return v0;
            yield return v1;
            yield return v2;
            yield return v3;
        }
        IEnumerator IEnumerable.GetEnumerator() { return GetEnumerator(); }
    }
    #endregion

    #region New class
    public static partial class New
    {
        public static T Type<T>() where T : new() { return new T(); }
        public static TResult Func<TResult>(Func<TResult> func) { return func(); }
        public static KeyValuePair<TKey, TValue> KeyValuePair<TKey, TValue>(TKey key, TValue value) { return new KeyValuePair<TKey, TValue>(key, value); }
        public static V2<T0, T1> V2<T0, T1>(T0 v0, T1 v1) { return new V2<T0, T1>(v0, v1); }
        public static V3<T0, T1, T2> V3<T0, T1, T2>(T0 v0, T1 v1, T2 v2) { return new V3<T0, T1, T2>(v0, v1, v2); }
        public static V4<T0, T1, T2, T3> V4<T0, T1, T2, T3>(T0 v0, T1 v1, T2 v2, T3 v3) { return new V4<T0, T1, T2, T3>(v0, v1, v2, v3); }
        public static Array1<T> Array1<T>(T v0) { return new Array1<T>(v0); }
        public static Array2<T> Array2<T>(T v0, T v1) { return new Array2<T>(v0, v1); }

        public static TResult[] Array<TResult>(int count, Func<int, TResult> func)
        {
            TResult[] array = new TResult[count];
            for (int i = 0; i < count; i++)
                array[i] = func(i);
            return array;
        }
        public static TResult[,] Array<TResult>(int count0, int count1, Func<int, int, TResult> func)
        {
            TResult[,] array = new TResult[count0, count1];
            for (int i0 = 0; i0 < count0; i0++)
                for (int i1 = 0; i1 < count1; i1++)
                    array[i0, i1] = func(i0, i1);
            return array;
        }
        public static TResult[, ,] Array<TResult>(int count0, int count1, int count2, Func<int, int, int, TResult> func)
        {
            TResult[, ,] array = new TResult[count0, count1, count2];
            for (int i0 = 0; i0 < count0; i0++)
                for (int i1 = 0; i1 < count1; i1++)
                    for (int i2 = 0; i2 < count2; i2++)
                        array[i0, i1, i2] = func(i0, i1, i2);
            return array;
        }
        public static TResult[, , ,] Array<TResult>(int count0, int count1, int count2, int count3, Func<int, int, int, int, TResult> func)
        {
            TResult[, , ,] array = new TResult[count0, count1, count2, count3];
            for (int i0 = 0; i0 < count0; i0++)
                for (int i1 = 0; i1 < count1; i1++)
                    for (int i2 = 0; i2 < count2; i2++)
                        for (int i3 = 0; i3 < count3; i3++)
                            array[i0, i1, i2, i3] = func(i0, i1, i2, i3);
            return array;
        }
        public static List<TResult> List<TResult>(int count, Func<int, TResult> func)
        {
            var list = new List<TResult>(count);
            for (int i = 0; i < count; i++)
                list.Add(func(i));
            return list;
        }

        public static TResult[] Array<TResult>(int count, Func<TResult> func)
        {
            TResult[] array = new TResult[count];
            for (int i = 0; i < count; i++)
                array[i] = func();
            return array;
        }
        public static TResult[,] Array<TResult>(int count0, int count1, Func<TResult> func)
        {
            TResult[,] array = new TResult[count0, count1];
            for (int i0 = 0; i0 < count0; i0++)
                for (int i1 = 0; i1 < count1; i1++)
                    array[i0, i1] = func();
            return array;
        }
        public static TResult[, ,] Array<TResult>(int count0, int count1, int count2, Func<TResult> func)
        {
            TResult[, ,] array = new TResult[count0, count1, count2];
            for (int i0 = 0; i0 < count0; i0++)
                for (int i1 = 0; i1 < count1; i1++)
                    for (int i2 = 0; i2 < count2; i2++)
                        array[i0, i1, i2] = func();
            return array;
        }
        public static TResult[, , ,] Array<TResult>(int count0, int count1, int count2, int count3, Func<TResult> func)
        {
            TResult[, , ,] array = new TResult[count0, count1, count2, count3];
            for (int i0 = 0; i0 < count0; i0++)
                for (int i1 = 0; i1 < count1; i1++)
                    for (int i2 = 0; i2 < count2; i2++)
                        for (int i3 = 0; i3 < count3; i3++)
                            array[i0, i1, i2, i3] = func();
            return array;
        }
        public static List<TResult> List<TResult>(int count, Func<TResult> func)
        {
            var list = new List<TResult>(count);
            for (int i = 0; i < count; i++)
                list.Add(func());
            return list;
        }

        public static TResult[] Array<TResult>(int count, TResult item)
        {
            TResult[] array = new TResult[count];
            for (int i = 0; i < count; i++)
                array[i] = item;
            return array;
        }
        public static TResult[,] Array<TResult>(int count0, int count1, TResult item)
        {
            TResult[,] array = new TResult[count0, count1];
            for (int i0 = 0; i0 < count0; i0++)
                for (int i1 = 0; i1 < count1; i1++)
                    array[i0, i1] = item;
            return array;
        }
        public static TResult[, ,] Array<TResult>(int count0, int count1, int count2, TResult item)
        {
            TResult[, ,] array = new TResult[count0, count1, count2];
            for (int i0 = 0; i0 < count0; i0++)
                for (int i1 = 0; i1 < count1; i1++)
                    for (int i2 = 0; i2 < count2; i2++)
                        array[i0, i1, i2] = item;
            return array;
        }
        public static TResult[, , ,] Array<TResult>(int count0, int count1, int count2, int count3, TResult item)
        {
            TResult[, , ,] array = new TResult[count0, count1, count2, count3];
            for (int i0 = 0; i0 < count0; i0++)
                for (int i1 = 0; i1 < count1; i1++)
                    for (int i2 = 0; i2 < count2; i2++)
                        for (int i3 = 0; i3 < count3; i3++)
                            array[i0, i1, i2, i3] = item;
            return array;
        }
        public static List<TResult> List<TResult>(int count, TResult item)
        {
            var list = new List<TResult>(count);
            for (int i = 0; i < count; i++)
                list.Add(item);
            return list;
        }

        public static StreamReader StreamReader(string path)
        {
            var encoding = Gb.DetectEncoding(path);
            if (encoding == null) return null;
            return new StreamReader(path, encoding);
        }
    }
    #endregion

    #region Int2
    public partial struct Int2 : IEquatable<Int2>, IComparable<Int2>
    {
        readonly int x;
        readonly int y;

        public Int2(int valueX, int valueY) { this.x = valueX; this.y = valueY; }
        public Int2(Int2 value) { this.x = value.x; this.y = value.y; }

        public int X { get { return x; } }
        public int Y { get { return y; } }

        public int this[int index]
        {
            get
            {
                switch (index)
                {
                    case 0: return x;
                    case 1: return y;
                    default: ThrowException.ArgumentOutOfRangeException("index"); return default(int);
                }
            }
        }

        public static bool operator ==(Int2 left, Int2 right) { return ((left.x - right.x) | (left.y - right.y)) == 0; }
        public static bool operator !=(Int2 left, Int2 right) { return ((left.x - right.x) | (left.y - right.y)) != 0; }
        public static bool operator <=(Int2 left, Int2 right) { return left.CompareTo(right) <= 0; }
        public static bool operator >=(Int2 left, Int2 right) { return left.CompareTo(right) >= 0; }
        public static bool operator <(Int2 left, Int2 right) { return left.CompareTo(right) < 0; }
        public static bool operator >(Int2 left, Int2 right) { return left.CompareTo(right) > 0; }
        public static Int2 operator ~(Int2 value) { return new Int2(~value.x, ~value.y); }
        public static Int2 operator -(Int2 value) { return new Int2(-value.x, -value.y); }
        public static Int2 operator +(Int2 left, Int2 right) { return new Int2(left.x + right.x, left.y + right.y); }
        public static Int2 operator -(Int2 left, Int2 right) { return new Int2(left.x - right.x, left.y - right.y); }
        public static Int2 operator *(Int2 left, Int2 right) { return new Int2(left.x * right.x, left.y * right.y); }
        public static Int2 operator /(Int2 left, Int2 right) { return new Int2(left.x / right.x, left.y / right.y); }
        public static Int2 operator %(Int2 left, Int2 right) { return new Int2(left.x % right.x, left.y % right.y); }
        public static Int2 operator |(Int2 left, Int2 right) { return new Int2(left.x | right.x, left.y | right.y); }
        public static Int2 operator &(Int2 left, Int2 right) { return new Int2(left.x & right.x, left.y & right.y); }
        public static Int2 operator ^(Int2 left, Int2 right) { return new Int2(left.x ^ right.x, left.y ^ right.y); }
        public static Int2 operator |(Int2 left, int right) { return new Int2(left.x | right, left.y | right); }
        public static Int2 operator &(Int2 left, int right) { return new Int2(left.x & right, left.y & right); }
        public static Int2 operator ^(Int2 left, int right) { return new Int2(left.x ^ right, left.y ^ right); }
        public static Int2 operator +(Int2 left, int right) { return new Int2(left.x + right, left.y + right); }
        public static Int2 operator -(Int2 left, int right) { return new Int2(left.x - right, left.y - right); }
        public static Int2 operator *(Int2 left, int right) { return new Int2(left.x * right, left.y * right); }
        public static Int2 operator /(Int2 left, int right) { return new Int2(left.x / right, left.y / right); }
        public static Int2 operator %(Int2 left, int right) { return new Int2(left.x % right, left.y % right); }

        public static Int2 OnesComplement(Int2 value) { return new Int2(~value.x, ~value.y); }
        public static Int2 Negate(Int2 value) { return new Int2(-value.x, -value.y); }
        public static Int2 Add(Int2 left, Int2 right) { return new Int2(left.x + right.x, left.y + right.y); }
        public static Int2 Subtract(Int2 left, Int2 right) { return new Int2(left.x - right.x, left.y - right.y); }
        public static Int2 Multiply(Int2 left, Int2 right) { return new Int2(left.x * right.x, left.y * right.y); }
        public static Int2 Divide(Int2 left, Int2 right) { return new Int2(left.x / right.x, left.y / right.y); }
        public static Int2 Mod(Int2 left, Int2 right) { return new Int2(left.x % right.x, left.y % right.y); }
        public static Int2 BitwiseOr(Int2 left, Int2 right) { return new Int2(left.x | right.x, left.y | right.y); }
        public static Int2 BitwiseAnd(Int2 left, Int2 right) { return new Int2(left.x & right.x, left.y & right.y); }
        public static Int2 Xor(Int2 left, Int2 right) { return new Int2(left.x ^ right.x, left.y ^ right.y); }
        public static Int2 BitwiseOr(Int2 left, int right) { return new Int2(left.x | right, left.y | right); }
        public static Int2 BitwiseAnd(Int2 left, int right) { return new Int2(left.x & right, left.y & right); }
        public static Int2 Xor(Int2 left, int right) { return new Int2(left.x ^ right, left.y ^ right); }
        public static Int2 Add(Int2 left, int right) { return new Int2(left.x + right, left.y + right); }
        public static Int2 Subtract(Int2 left, int right) { return new Int2(left.x - right, left.y - right); }
        public static Int2 Multiply(Int2 left, int right) { return new Int2(left.x * right, left.y * right); }
        public static Int2 Divide(Int2 left, int right) { return new Int2(left.x / right, left.y / right); }
        public static Int2 Mod(Int2 left, int right) { return new Int2(left.x % right, left.y % right); }

        public override bool Equals(object obj) { return (obj is Int2) && (this == (Int2)obj); }
        public override int GetHashCode() { return ((y << 16) | (int)((uint)y >> 16)) ^ x; }
        public override string ToString() { return string.Format("{0}, {1}", x, y); }

        #region IEquatable<Int2> メンバ
        public bool Equals(Int2 other) { return ((x - other.x) | (y - other.y)) == 0; }
        #endregion
        #region IComparable<Int2> メンバ
        public int CompareTo(Int2 other) { return y != other.y ? y - other.y : x - other.x; }
        #endregion
    }
    #endregion

    #region Int3
    public partial struct Int3 : IEquatable<Int3>, IComparable<Int3>
    {
        readonly int x;
        readonly int y;
        readonly int z;

        public Int3(int valueX, int valueY, int valueZ) { this.x = valueX; this.y = valueY; this.z = valueZ; }
        public Int3(Int3 value) { this.x = value.x; this.y = value.y; this.z = value.z; }

        public int X { get { return x; } }
        public int Y { get { return y; } }
        public int Z { get { return z; } }

        public int this[int index]
        {
            get
            {
                switch (index)
                {
                    case 0: return x;
                    case 1: return y;
                    case 2: return z;
                    default: ThrowException.ArgumentOutOfRangeException("index"); return default(int);
                }
            }
        }

        public static bool operator ==(Int3 left, Int3 right) { return ((left.x - right.x) | (left.y - right.y) | (left.z - right.z)) == 0; }
        public static bool operator !=(Int3 left, Int3 right) { return ((left.x - right.x) | (left.y - right.y) | (left.z - right.z)) != 0; }
        public static bool operator <=(Int3 left, Int3 right) { return left.CompareTo(right) <= 0; }
        public static bool operator >=(Int3 left, Int3 right) { return left.CompareTo(right) >= 0; }
        public static bool operator <(Int3 left, Int3 right) { return left.CompareTo(right) < 0; }
        public static bool operator >(Int3 left, Int3 right) { return left.CompareTo(right) > 0; }
        public static Int3 operator ~(Int3 value) { return new Int3(~value.x, ~value.y, ~value.z); }
        public static Int3 operator -(Int3 value) { return new Int3(-value.x, -value.y, -value.z); }
        public static Int3 operator +(Int3 left, Int3 right) { return new Int3(left.x + right.x, left.y + right.y, left.z + right.z); }
        public static Int3 operator -(Int3 left, Int3 right) { return new Int3(left.x - right.x, left.y - right.y, left.z - right.z); }
        public static Int3 operator *(Int3 left, Int3 right) { return new Int3(left.x * right.x, left.y * right.y, left.z * right.z); }
        public static Int3 operator /(Int3 left, Int3 right) { return new Int3(left.x / right.x, left.y / right.y, left.z / right.z); }
        public static Int3 operator %(Int3 left, Int3 right) { return new Int3(left.x % right.x, left.y % right.y, left.z % right.z); }
        public static Int3 operator |(Int3 left, Int3 right) { return new Int3(left.x | right.x, left.y | right.y, left.z | right.z); }
        public static Int3 operator &(Int3 left, Int3 right) { return new Int3(left.x & right.x, left.y & right.y, left.z & right.z); }
        public static Int3 operator ^(Int3 left, Int3 right) { return new Int3(left.x ^ right.x, left.y ^ right.y, left.z ^ right.z); }
        public static Int3 operator |(Int3 left, int right) { return new Int3(left.x | right, left.y | right, left.z | right); }
        public static Int3 operator &(Int3 left, int right) { return new Int3(left.x & right, left.y & right, left.z & right); }
        public static Int3 operator ^(Int3 left, int right) { return new Int3(left.x ^ right, left.y ^ right, left.z ^ right); }
        public static Int3 operator +(Int3 left, int right) { return new Int3(left.x + right, left.y + right, left.z + right); }
        public static Int3 operator -(Int3 left, int right) { return new Int3(left.x - right, left.y - right, left.z - right); }
        public static Int3 operator *(Int3 left, int right) { return new Int3(left.x * right, left.y * right, left.z * right); }
        public static Int3 operator /(Int3 left, int right) { return new Int3(left.x / right, left.y / right, left.z / right); }
        public static Int3 operator %(Int3 left, int right) { return new Int3(left.x % right, left.y % right, left.z % right); }

        public static Int3 OnesComplement(Int3 value) { return new Int3(~value.x, ~value.y, ~value.z); }
        public static Int3 Negate(Int3 value) { return new Int3(-value.x, -value.y, -value.z); }
        public static Int3 Add(Int3 left, Int3 right) { return new Int3(left.x + right.x, left.y + right.y, left.z + right.z); }
        public static Int3 Subtract(Int3 left, Int3 right) { return new Int3(left.x - right.x, left.y - right.y, left.z - right.z); }
        public static Int3 Multiply(Int3 left, Int3 right) { return new Int3(left.x * right.x, left.y * right.y, left.z * right.z); }
        public static Int3 Divide(Int3 left, Int3 right) { return new Int3(left.x / right.x, left.y / right.y, left.z / right.z); }
        public static Int3 Mod(Int3 left, Int3 right) { return new Int3(left.x % right.x, left.y % right.y, left.z % right.z); }
        public static Int3 BitwiseOr(Int3 left, Int3 right) { return new Int3(left.x | right.x, left.y | right.y, left.z | right.z); }
        public static Int3 BitwiseAnd(Int3 left, Int3 right) { return new Int3(left.x & right.x, left.y & right.y, left.z & right.z); }
        public static Int3 Xor(Int3 left, Int3 right) { return new Int3(left.x ^ right.x, left.y ^ right.y, left.z ^ right.z); }
        public static Int3 BitwiseOr(Int3 left, int right) { return new Int3(left.x | right, left.y | right, left.z | right); }
        public static Int3 BitwiseAnd(Int3 left, int right) { return new Int3(left.x & right, left.y & right, left.z & right); }
        public static Int3 Xor(Int3 left, int right) { return new Int3(left.x ^ right, left.y ^ right, left.z ^ right); }
        public static Int3 Add(Int3 left, int right) { return new Int3(left.x + right, left.y + right, left.z + right); }
        public static Int3 Subtract(Int3 left, int right) { return new Int3(left.x - right, left.y - right, left.z - right); }
        public static Int3 Multiply(Int3 left, int right) { return new Int3(left.x * right, left.y * right, left.z * right); }
        public static Int3 Divide(Int3 left, int right) { return new Int3(left.x / right, left.y / right, left.z / right); }
        public static Int3 Mod(Int3 left, int right) { return new Int3(left.x % right, left.y % right, left.z % right); }

        public override bool Equals(object obj) { return (obj is Int3) && (this == (Int3)obj); }
        public override int GetHashCode() { return ((y << 16) | (int)((uint)y >> 16)) ^ x ^ z; }
        public override string ToString() { return string.Format("{0}, {1}, {2}", x, y, z); }

        #region IEquatable<Int3> メンバ
        public bool Equals(Int3 other) { return ((x - other.x) | (y - other.y) | (z - other.z)) == 0; }
        #endregion
        #region IComparable<Int3> メンバ
        public int CompareTo(Int3 other) { return z != other.z ? z - other.z : (y != other.y ? y - other.y : x - other.x); }
        #endregion
    }
    #endregion

    #region Double2
    public partial struct Double2 : IEquatable<Double2>, IComparable<Double2>
    {
        readonly double x;
        readonly double y;

        public Double2(double valueX, double valueY) { this.x = valueX; this.y = valueY; }
        public Double2(Double2 value) { this.x = value.x; this.y = value.y; }
        public Double2(Int2 value) { this.x = value.X; this.y = value.Y; }

        public double X { get { return x; } }
        public double Y { get { return y; } }

        public double this[int index]
        {
            get
            {
                switch (index)
                {
                    case 0: return x;
                    case 1: return y;
                    default: ThrowException.ArgumentOutOfRangeException("index"); return default(double);
                }
            }
        }

        public static bool operator ==(Double2 left, Double2 right) { return (left.x == right.x) && (left.y == right.y); }
        public static bool operator !=(Double2 left, Double2 right) { return (left.x != right.x) || (left.y != right.y); }
        public static bool operator <=(Double2 left, Double2 right) { return left.CompareTo(right) <= 0; }
        public static bool operator >=(Double2 left, Double2 right) { return left.CompareTo(right) >= 0; }
        public static bool operator <(Double2 left, Double2 right) { return left.CompareTo(right) < 0; }
        public static bool operator >(Double2 left, Double2 right) { return left.CompareTo(right) > 0; }
        public static Double2 operator -(Double2 value) { return new Double2(-value.x, -value.y); }
        public static Double2 operator +(Double2 left, Double2 right) { return new Double2(left.x + right.x, left.y + right.y); }
        public static Double2 operator -(Double2 left, Double2 right) { return new Double2(left.x - right.x, left.y - right.y); }
        public static Double2 operator *(Double2 left, Double2 right) { return new Double2(left.x * right.x, left.y * right.y); }
        public static Double2 operator /(Double2 left, Double2 right) { return new Double2(left.x / right.x, left.y / right.y); }
        public static Double2 operator +(Double2 left, double right) { return new Double2(left.x + right, left.y + right); }
        public static Double2 operator -(Double2 left, double right) { return new Double2(left.x - right, left.y - right); }
        public static Double2 operator *(Double2 left, double right) { return new Double2(left.x * right, left.y * right); }
        public static Double2 operator /(Double2 left, double right) { return new Double2(left.x / right, left.y / right); }

        public static Double2 Negate(Double2 value) { return new Double2(-value.x, -value.y); }
        public static Double2 Add(Double2 left, Double2 right) { return new Double2(left.x + right.x, left.y + right.y); }
        public static Double2 Subtract(Double2 left, Double2 right) { return new Double2(left.x - right.x, left.y - right.y); }
        public static Double2 Multiply(Double2 left, Double2 right) { return new Double2(left.x * right.x, left.y * right.y); }
        public static Double2 Divide(Double2 left, Double2 right) { return new Double2(left.x / right.x, left.y / right.y); }
        public static Double2 Add(Double2 left, double right) { return new Double2(left.x + right, left.y + right); }
        public static Double2 Subtract(Double2 left, double right) { return new Double2(left.x - right, left.y - right); }
        public static Double2 Multiply(Double2 left, double right) { return new Double2(left.x * right, left.y * right); }
        public static Double2 Divide(Double2 left, double right) { return new Double2(left.x / right, left.y / right); }

        public static bool IsNaN(Double2 value) { return double.IsNaN(value.x) || double.IsNaN(value.y); }
        public static readonly Double2 Zero = new Double2(0, 0);
        public static readonly Double2 One = new Double2(1, 1);

        public override bool Equals(object obj) { return (obj is Double2) && (this == (Double2)obj); }
        public override int GetHashCode() { return (int)(y - x); }
        public override string ToString() { return string.Format("{0}, {1}", x, y); }

        #region IEquatable<Double2> メンバ
        public bool Equals(Double2 other) { return (x == other.x) && (y == other.y); }
        #endregion
        #region IComparable<Double2> メンバ
        public int CompareTo(Double2 other) { return y != other.y ? (y < other.y ? -1 : 1) : (x == other.x ? 0 : x < other.x ? -1 : 1); }
        #endregion
    }
    #endregion

    #region Double3
    public partial struct Double3 : IEquatable<Double3>, IComparable<Double3>
    {
        readonly double x;
        readonly double y;
        readonly double z;

        public Double3(double valueX, double valueY, double valueZ) { this.x = valueX; this.y = valueY; this.z = valueZ; }
        public Double3(Double3 value) { this.x = value.x; this.y = value.y; this.z = value.z; }
        public Double3(Int3 value) { this.x = value.X; this.y = value.Y; this.z = value.Z; }

        public double X { get { return x; } }
        public double Y { get { return y; } }
        public double Z { get { return z; } }

        public double this[int index]
        {
            get
            {
                switch (index)
                {
                    case 0: return x;
                    case 1: return y;
                    case 2: return z;
                    default: ThrowException.ArgumentOutOfRangeException("index"); return default(double);
                }
            }
        }
        public static bool operator ==(Double3 left, Double3 right) { return (left.x == right.x) && (left.y == right.y) && (left.z == right.z); }
        public static bool operator !=(Double3 left, Double3 right) { return (left.x != right.x) || (left.y != right.y) || (left.z != right.z); }
        public static bool operator <=(Double3 left, Double3 right) { return left.CompareTo(right) <= 0; }
        public static bool operator >=(Double3 left, Double3 right) { return left.CompareTo(right) >= 0; }
        public static bool operator <(Double3 left, Double3 right) { return left.CompareTo(right) < 0; }
        public static bool operator >(Double3 left, Double3 right) { return left.CompareTo(right) > 0; }
        public static Double3 operator -(Double3 value) { return new Double3(-value.x, -value.y, -value.z); }
        public static Double3 operator +(Double3 left, Double3 right) { return new Double3(left.x + right.x, left.y + right.y, left.z + right.z); }
        public static Double3 operator -(Double3 left, Double3 right) { return new Double3(left.x - right.x, left.y - right.y, left.z - right.z); }
        public static Double3 operator *(Double3 left, Double3 right) { return new Double3(left.x * right.x, left.y * right.y, left.z * right.z); }
        public static Double3 operator /(Double3 left, Double3 right) { return new Double3(left.x / right.x, left.y / right.y, left.z / right.z); }
        public static Double3 operator +(Double3 left, double right) { return new Double3(left.x + right, left.y + right, left.z + right); }
        public static Double3 operator -(Double3 left, double right) { return new Double3(left.x - right, left.y - right, left.z - right); }
        public static Double3 operator *(Double3 left, double right) { return new Double3(left.x * right, left.y * right, left.z * right); }
        public static Double3 operator /(Double3 left, double right) { return new Double3(left.x / right, left.y / right, left.z / right); }

        public static Double3 Negate(Double3 value) { return new Double3(-value.x, -value.y, -value.z); }
        public static Double3 Add(Double3 left, Double3 right) { return new Double3(left.x + right.x, left.y + right.y, left.z + right.z); }
        public static Double3 Subtract(Double3 left, Double3 right) { return new Double3(left.x - right.x, left.y - right.y, left.z - right.z); }
        public static Double3 Multiply(Double3 left, Double3 right) { return new Double3(left.x * right.x, left.y * right.y, left.z * right.z); }
        public static Double3 Divide(Double3 left, Double3 right) { return new Double3(left.x / right.x, left.y / right.y, left.z / right.z); }
        public static Double3 Add(Double3 left, double right) { return new Double3(left.x + right, left.y + right, left.z + right); }
        public static Double3 Subtract(Double3 left, double right) { return new Double3(left.x - right, left.y - right, left.z - right); }
        public static Double3 Multiply(Double3 left, double right) { return new Double3(left.x * right, left.y * right, left.z * right); }
        public static Double3 Divide(Double3 left, double right) { return new Double3(left.x / right, left.y / right, left.z / right); }

        public static bool IsNaN(Double3 value) { return double.IsNaN(value.x) || double.IsNaN(value.y) || double.IsNaN(value.z); }
        public static readonly Double3 Zero = new Double3(0, 0, 0);
        public static readonly Double3 One = new Double3(1, 1, 1);

        public override bool Equals(object obj) { return (obj is Double3) && (this == (Double3)obj); }
        public override int GetHashCode() { return (int)(y - x + z); }
        public override string ToString() { return string.Format("{0}, {1}, {2}", x, y, z); }

        #region IEquatable<Double3> メンバ
        public bool Equals(Double3 other) { return (x == other.x) && (y == other.y) && (z == other.z); }
        #endregion
        #region IComparable<Double3> メンバ
        public int CompareTo(Double3 other) { return z != other.z ? (z < other.y ? -1 : 1) : (y != other.y ? (y < other.y ? -1 : 1) : (x == other.x ? 0 : x < other.x ? -1 : 1)); }
        #endregion
    }
    #endregion

    #region Ints
    class Ints : IComparable, IEquatable<Ints>, IComparable<Ints>, IList<int>
    {
        readonly int[] _Items;
        public Ints(params int[] items)
        {
            this._Items = items.CloneX();
        }
        public Ints(IEnumerable<int> collections)
        {
            this._Items = collections.ToArray();
        }

        public int Length
        {
            get { return _Items.Length; }
        }

        public static bool operator ==(Ints left, Ints right) { return Equals(left, right); }
        public static bool operator !=(Ints left, Ints right) { return !Equals(left, right); }
        public static bool operator <=(Ints left, Ints right) { return Compare(left, right) <= 0; }
        public static bool operator >=(Ints left, Ints right) { return Compare(left, right) >= 0; }
        public static bool operator <(Ints left, Ints right) { return Compare(left, right) < 0; }
        public static bool operator >(Ints left, Ints right) { return Compare(left, right) > 0; }
        public static Ints operator +(Ints x, Ints y) { return new Ints(x._Items.Concat(y._Items)); }

        public static int Compare(Ints left, Ints right)
        {
            if ((object)left == null) return ((object)right == null) ? 0 : -1;
            return left.CompareTo(right);
        }
        public static bool Equals(Ints left, Ints right)
        {
            if ((object)left == null) return (object)right == null;
            return left.Equals(right);
        }

        public override bool Equals(object obj) { return (obj is Ints) && Equals((Ints)obj); }
        public override int GetHashCode()
        {
            int h = 0;
            for (int i = 0; i < _Items.Length; i++)
                h = h * 3 - 13 + _Items[i];
            return h;
        }
        public override string ToString() { return string.Join(",", _Items.Select(l => l.ToString()).ToArray()); }

        #region IComparable メンバー
        public int CompareTo(object obj)
        {
            var other = obj as Ints;
            if (other == null) ThrowException.ArgumentException("obj");
            return CompareTo(other);
        }
        #endregion
        #region IEquatable<Ints> メンバ
        public bool Equals(Ints other)
        {
            if ((object)other == null) return false;
            if (_Items.Length != other._Items.Length) return false;
            for (int i = 0; i < _Items.Length; i++)
                if (_Items[i] != other._Items[i]) return false;
            return true;
        }
        #endregion
        #region IComparable<Ints> メンバ
        public int CompareTo(Ints other)
        {
            if ((object)other == null) return 1;
            int l = Math.Min(_Items.Length, other._Items.Length);
            for (int i = 0; i < l; i++)
                if (_Items[i] != other._Items[i]) return _Items[i] - other._Items[i];
            return _Items.Length - other._Items.Length;
        }
        #endregion

        #region IList<int> メンバー
        public int IndexOf(int item) { return _Items.IndexOf(item); }
        public void Insert(int index, int item) { throw new NotImplementedException(); }
        public void RemoveAt(int index) { throw new NotImplementedException(); }
        public int this[int index]
        {
            get { return _Items[index]; }
            set { throw new NotImplementedException(); }
        }
        #endregion
        #region ICollection<int> メンバー
        public void Add(int item) { throw new NotImplementedException(); }
        public void Clear() { throw new NotImplementedException(); }
        public bool Contains(int item) { return _Items.Contains(item); }
        public void CopyTo(int[] array, int arrayIndex) { _Items.CopyTo(array, arrayIndex); }
        public int Count
        {
            get { return _Items.Length; }
        }
        public bool IsReadOnly
        {
            get { return true; }
        }
        public bool Remove(int item) { throw new NotImplementedException(); }
        #endregion
        #region IEnumerable<int> メンバー
        public IEnumerator<int> GetEnumerator()
        {
            for (int i = 0; i < _Items.Length; i++)
                yield return _Items[i];
        }
        #endregion
        #region IEnumerable メンバー
        IEnumerator IEnumerable.GetEnumerator() { return GetEnumerator(); }
        #endregion
    }
    #endregion

    #region Extension class
    public static partial class Ex
    {
        #region code snippet
        public static IEnumerable<int> Range(int count) { return Enumerable.Range(0, count); }
        public static IEnumerable<int> FromTo(int start, int end)
        {
            if (start <= end)
                for (int i = start; i < end; i++)
                    yield return i;
            else
                for (int i = start; --i >= end; )
                    yield return i;
        }
        public static IEnumerable<int> FromUntil(int start, int end)
        {
            if (start <= end)
                for (int i = start; i <= end; i++)
                    yield return i;
            else
                for (int i = start; i >= end; i--)
                    yield return i;
        }
        #endregion

        #region Array
        public static int IndexOf<T>(this T[] array, T value) { return Array.IndexOf(array, value); }
        public static void Clear(this Array array)
        {
            Array.Clear(array, 0, array.Length);
        }

        public static void CopyTo(this Array sourceArray, Array destinationArray)
        {
            sourceArray.CopyTo(destinationArray, 0);
        }
        public static void CopyTo(this Array sourceArray, int srcIndex, Array destinationArray, int dstinationIndex, int length)
        {
            Array.Copy(sourceArray, srcIndex, destinationArray, dstinationIndex, length);
        }
        public static T[] CloneX<T>(this T[] array)
        {
            return (T[])array.Clone();
        }
        public static T[,] CloneX<T>(this T[,] array)
        {
            return (T[,])array.Clone();
        }
        public static T[, ,] CloneX<T>(this T[, ,] array)
        {
            return (T[, ,])array.Clone();
        }
        public static T[, , ,] CloneX<T>(this T[, , ,] array)
        {
            return (T[, , ,])array.Clone();
        }
        public static T[][] DeepClone<T>(this T[][] array)
        {
            return New.Array(array.Length, i => array[i].CloneX());
        }
        public static T[][][] DeepClone<T>(this T[][][] array)
        {
            return New.Array(array.Length, i => array[i].DeepClone());
        }
        public static T[][][][] DeepClone<T>(this T[][][][] array)
        {
            return New.Array(array.Length, i => array[i].DeepClone());
        }

        public static void Let<T>(this T[] array, T value)
        {
            for (int i = array.Length; --i >= 0; )
                array[i] = value;
        }
        public static void Let<T>(this T[,] array, T value)
        {
            for (int i = array.GetLength(0); --i >= 0; )
                for (int j = array.GetLength(1); --j >= 0; )
                    array[i, j] = value;
        }
        public static void Let<T>(this T[, ,] array, T value)
        {
            for (int i = array.GetLength(0); --i >= 0; )
                for (int j = array.GetLength(1); --j >= 0; )
                    for (int k = array.GetLength(2); --k >= 0; )
                        array[i, j, k] = value;
        }
        public static void Let<T>(this T[, , ,] array, T value)
        {
            for (int i = array.GetLength(0); --i >= 0; )
                for (int j = array.GetLength(1); --j >= 0; )
                    for (int k = array.GetLength(2); --k >= 0; )
                        for (int l = array.GetLength(3); --l >= 0; )
                            array[i, j, k, l] = value;
        }
        public static T[] Sub<T>(this T[] array, int start)
        {
            T[] data = new T[array.Length - start];
            Array.Copy(array, start, data, 0, data.Length);
            return data;
        }
        public static T[] Sub<T>(this T[] array, int start, int count)
        {
            T[] data = new T[count];
            Array.Copy(array, start, data, 0, count);
            return data;
        }

        public static T[][] ToJaggedArray<T>(this T[,] array)
        {
            var jagged = new T[array.GetLength(0)][];
            for (int i = 0; i < jagged.Length; i++)
            {
                var row = new T[array.GetLength(1)];
                for (int j = 0; j < row.Length; j++)
                    row[j] = array[i, j];
                jagged[i] = row;
            }
            return jagged;
        }
        #endregion

        #region List
        public static bool AddOrDiscard<T>(this List<T> list, T item)
        {
            if (list.Contains(item)) return false;
            list.Add(item); return true;
        }
        public static bool AddOrOverwrite<T>(this List<T> list, T item)
        {
            int index = list.IndexOf(item);
            if (index != -1) { list[index] = item; return false; }
            list.Add(item); return true;
        }
        #endregion

        #region IEnumerable
        public static TSource FirstOrDefault<TSource>(this IEnumerable<TSource> source, TSource value, Func<TSource, bool> predicate)
        {
            foreach (var element in source)
                if (predicate(element)) return element;
            return value;
        }
        public static void ForEach<TSource>(this IEnumerable<TSource> source, Action<TSource> action)
        {
            foreach (var element in source) action(element);
        }
        public static void ForEach<TSource>(this IEnumerable<TSource> source, Action<TSource, int> action)
        {
            int count = 0;
            foreach (var element in source) action(element, count++);
        }
        public static IOrderedEnumerable<TSource> OrderBy<TSource>(this IEnumerable<TSource> source) { return source.OrderBy(x => x); }
        public static IOrderedEnumerable<TSource> OrderByDescending<TSource>(this IEnumerable<TSource> source) { return source.OrderByDescending(x => x); }

        public static Dictionary<TKey, TValue> ToDictionary<TKey, TValue>(this IEnumerable<KeyValuePair<TKey, TValue>> source)
        {
            var dictionary = new Dictionary<TKey, TValue>();
            foreach (var pair in source)
                dictionary.Add(pair.Key, pair.Value);
            return dictionary;
        }
        public static IEnumerable<T> Row<T>(this T[,] matrix, int row)
        {
            return Enumerable.Range(0, matrix.GetLength(1)).Select(j => matrix[row, j]);
        }
        public static IEnumerable<T> Col<T>(this T[,] matrix, int col)
        {
            return Enumerable.Range(0, matrix.GetLength(0)).Select(i => matrix[i, col]);
        }

        public static Dictionary<TKey, int> ToHistogram<TKey>(this IEnumerable<TKey> source)
        {
            var dictionary = new Dictionary<TKey, int>();
            foreach (var element in source)
                if (dictionary.ContainsKey(element)) dictionary[element]++;
                else dictionary.Add(element, 1);
            return dictionary;
        }

        public static Dictionary<TKey, int> ToHistogram<TKey>(this IEnumerable<TKey> source, IEqualityComparer<TKey> comparer)
        {
            var dictionary = new Dictionary<TKey, int>(comparer);
            foreach (var element in source)
                if (dictionary.ContainsKey(element)) dictionary[element]++;
                else dictionary.Add(element, 1);
            return dictionary;
        }

        public static IEnumerable<string> SelectString<TSource>(this IEnumerable<TSource> source) { return source.Select(element => element.ToString()); }
        public static IEnumerable<string> SelectString(this IEnumerable<byte> source, string format) { return source.Select(element => element.ToString(format)); }
        public static IEnumerable<string> SelectString(this IEnumerable<sbyte> source, string format) { return source.Select(element => element.ToString(format)); }
        public static IEnumerable<string> SelectString(this IEnumerable<ushort> source, string format) { return source.Select(element => element.ToString(format)); }
        public static IEnumerable<string> SelectString(this IEnumerable<short> source, string format) { return source.Select(element => element.ToString(format)); }
        public static IEnumerable<string> SelectString(this IEnumerable<uint> source, string format) { return source.Select(element => element.ToString(format)); }
        public static IEnumerable<string> SelectString(this IEnumerable<int> source, string format) { return source.Select(element => element.ToString(format)); }
        public static IEnumerable<string> SelectString(this IEnumerable<ulong> source, string format) { return source.Select(element => element.ToString(format)); }
        public static IEnumerable<string> SelectString(this IEnumerable<long> source, string format) { return source.Select(element => element.ToString(format)); }
        public static IEnumerable<string> SelectString(this IEnumerable<float> source, string format) { return source.Select(element => element.ToString(format)); }
        public static IEnumerable<string> SelectString(this IEnumerable<double> source, string format) { return source.Select(element => element.ToString(format)); }
        public static IEnumerable<string> SelectString(this IEnumerable<decimal> source, string format) { return source.Select(element => element.ToString(format)); }
        public static IEnumerable<string> SelectString(this IEnumerable<DateTime> source, string format) { return source.Select(element => element.ToString(format)); }

        #region Calculation
        public static uint Sum<TSource>(this IEnumerable<TSource> source, Func<TSource, uint> selector) { return source.Select(element => selector(element)).Sum(); }
        public static uint Sum(this IEnumerable<uint> source)
        {
            uint sum = 0;
            foreach (var element in source) checked { sum += element; }
            return sum;
        }
        public static ulong Sum<TSource>(this IEnumerable<TSource> source, Func<TSource, ulong> selector) { return source.Select(element => selector(element)).Sum(); }
        public static ulong Sum(this IEnumerable<ulong> source)
        {
            ulong sum = 0;
            foreach (var element in source) checked { sum += element; }
            return sum;
        }

        public static double Product<TSource>(this IEnumerable<TSource> source, Func<TSource, double> selector) { return source.Select(element => selector(element)).Product(); }
        public static double Product(this IEnumerable<double> source)
        {
            double product = 1;
            foreach (var element in source) checked { product *= element; }
            return product;
        }

        public static double Variance<TSource>(this IEnumerable<TSource> source, Func<TSource, double> selector) { return source.Select(element => selector(element)).Variance(); }
        public static double Variance(this IEnumerable<double> source)
        {
            int count = 0;
            double moment1 = 0;
            double moment2 = 0;
            foreach (var element in source)
            {
                checked
                {
                    moment1 += element;
                    moment2 += element * element;
                    count++;
                }
            }
            if (count <= 1) ThrowException.InvalidOperationException("count <= 1");
            return Math.Max(0.0, moment2 - moment1 * moment1 / count) / (count - 1);
        }

        public static double VariancePopulation<TSource>(this IEnumerable<TSource> source, Func<TSource, double> selector) { return source.Select(element => selector(element)).VariancePopulation(); }
        public static double VariancePopulation(this IEnumerable<double> source)
        {
            int count = 0;
            double moment1 = 0;
            double moment2 = 0;
            foreach (var element in source)
            {
                checked
                {
                    moment1 += element;
                    moment2 += element * element;
                    count++;
                }
            }
            if (count <= 0) ThrowException.NoElements();
            return Math.Max(0.0, moment2 - moment1 * moment1 / count) / count;
        }

        public static int[] Cumurative(this IEnumerable<int> source)
        {
            int sum = 0;
            return source.Select(element => checked(sum += element)).ToArray();
        }
        public static double[] Cumurative(this IEnumerable<double> source)
        {
            double sum = 0;
            return source.Select(element => checked(sum += element)).ToArray();
        }
        public static int[] CumurativeZero(this IEnumerable<int> source)
        {
            int sum = 0;
            return source.Select(element => { var s = sum; checked { sum += element; } return s; }).ToArray();
        }
        public static double[] CumurativeZero(this IEnumerable<double> source)
        {
            double sum = 0;
            return source.Select(element => { var s = sum; checked { sum += element; } return s; }).ToArray();
        }
        #endregion

        #region MinMax
        public static V2<V2<int, TSource>, V2<int, TSource>> MinMaxIndexItem<TSource>(this IEnumerable<TSource> source)
            where TSource : IComparable<TSource>
        {
            var minIndex = default(int);
            var maxIndex = default(int);
            var minItem = default(TSource);
            var maxItem = default(TSource);
            int index = 0;
            foreach (var item in source)
            {
                if (index == 0 || minItem.CompareTo(item) > 0) { minItem = item; minIndex = index; }
                if (index == 0 || maxItem.CompareTo(item) < 0) { maxItem = item; maxIndex = index; }
                index++;
            }
            if (index == 0) ThrowException.NoElements();
            return New.V2(New.V2(minIndex, minItem), New.V2(maxIndex, maxItem));
        }
        public static V2<int, TSource> MinIndexItem<TSource>(this IEnumerable<TSource> source)
            where TSource : IComparable<TSource>
        {
            var minIndex = default(int);
            var minItem = default(TSource);
            int index = 0;
            foreach (var item in source)
            {
                if (index == 0 || minItem.CompareTo(item) > 0) { minItem = item; minIndex = index; }
                index++;
            }
            if (index == 0) ThrowException.NoElements();
            return New.V2(minIndex, minItem);
        }
        public static V2<int, TSource> MaxIndexItem<TSource>(this IEnumerable<TSource> source)
            where TSource : IComparable<TSource>
        {
            var maxIndex = default(int);
            var maxItem = default(TSource);
            int index = 0;
            foreach (var item in source)
            {
                if (index == 0 || maxItem.CompareTo(item) < 0) { maxItem = item; maxIndex = index; }
                index++;
            }
            if (index == 0) ThrowException.NoElements();
            return New.V2(maxIndex, maxItem);
        }
        public static int MinIndex<TSource>(this IEnumerable<TSource> source) where TSource : IComparable<TSource> { return source.MinIndexItem().v0; }
        public static int MaxIndex<TSource>(this IEnumerable<TSource> source) where TSource : IComparable<TSource> { return source.MaxIndexItem().v0; }

        public static V2<V3<int, TSource, TValue>, V3<int, TSource, TValue>> MinMaxIndexItemValue<TSource, TValue>(this IEnumerable<TSource> source, Func<TSource, TValue> selector)
            where TValue : IComparable<TValue>
        {
            var minIndex = default(int);
            var maxIndex = default(int);
            var minItem = default(TSource);
            var maxItem = default(TSource);
            var minValue = default(TValue);
            var maxValue = default(TValue);
            int index = 0;
            foreach (var item in source)
            {
                var value = selector(item);
                if (index == 0 || minValue.CompareTo(value) > 0) { minValue = value; minItem = item; minIndex = index; }
                if (index == 0 || maxValue.CompareTo(value) < 0) { maxValue = value; maxItem = item; maxIndex = index; }
                index++;
            }
            if (index == 0) ThrowException.NoElements();
            return New.V2(New.V3(minIndex, minItem, minValue), New.V3(maxIndex, maxItem, maxValue));
        }
        public static V3<int, TSource, TValue> MinIndexItemValue<TSource, TValue>(this IEnumerable<TSource> source, Func<TSource, TValue> selector)
            where TValue : IComparable<TValue>
        {
            var minIndex = default(int);
            var minItem = default(TSource);
            var minValue = default(TValue);
            int index = 0;
            foreach (var item in source)
            {
                var value = selector(item);
                if (index == 0 || minValue.CompareTo(value) > 0) { minValue = value; minItem = item; minIndex = index; }
                index++;
            }
            if (index == 0) ThrowException.NoElements();
            return New.V3(minIndex, minItem, minValue);
        }
        public static V3<int, TSource, TValue> MaxIndexItemValue<TSource, TValue>(this IEnumerable<TSource> source, Func<TSource, TValue> selector)
            where TValue : IComparable<TValue>
        {
            var maxIndex = default(int);
            var maxItem = default(TSource);
            var maxValue = default(TValue);
            int index = 0;
            foreach (var item in source)
            {
                var value = selector(item);
                if (index == 0 || maxValue.CompareTo(value) < 0) { maxValue = value; maxItem = item; maxIndex = index; }
                index++;
            }
            if (index == 0) ThrowException.NoElements();
            return New.V3(maxIndex, maxItem, maxValue);
        }
        public static int MinIndex<TSource, TValue>(this IEnumerable<TSource> source, Func<TSource, TValue> selector) where TValue : IComparable<TValue> { return source.MinIndexItemValue(selector).v0; }
        public static int MaxIndex<TSource, TValue>(this IEnumerable<TSource> source, Func<TSource, TValue> selector) where TValue : IComparable<TValue> { return source.MaxIndexItemValue(selector).v0; }
        public static TSource MinItem<TSource, TValue>(this IEnumerable<TSource> source, Func<TSource, TValue> selector) where TValue : IComparable<TValue> { return source.MinIndexItemValue(selector).v1; }
        public static TSource MaxItem<TSource, TValue>(this IEnumerable<TSource> source, Func<TSource, TValue> selector) where TValue : IComparable<TValue> { return source.MaxIndexItemValue(selector).v1; }
        #endregion
        #endregion

        #region IList
        // IList Index
        public static int MinIndex<T>(this IList<T> list, int index, int length) where T : IComparable<T>
        {
            int j;
            for (int i = j = index; ++i < index + length; )
                if (list[j].CompareTo(list[i]) > 0) j = i;
            return j;
        }
        public static int MaxIndex<T>(this IList<T> list, int index, int length) where T : IComparable<T>
        {
            int j;
            for (int i = j = index; ++i < index + length; )
                if (list[j].CompareTo(list[i]) < 0) j = i;
            return j;
        }
        public static Int2 MinMaxIndex<T>(this IList<T> list, int index, int length) where T : IComparable<T>
        {
            int j0, j1;
            for (int i = j0 = j1 = index; ++i < index + length; )
                if (list[j0].CompareTo(list[i]) > 0) j0 = i;
                else if (list[j1].CompareTo(list[i]) < 0) j1 = i;
            return new Int2(j0, j1);
        }
        public static int MinIndex<T>(this IList<T> list) where T : IComparable<T> { return MinIndex(list, 0, list.Count); }
        public static int MaxIndex<T>(this IList<T> list) where T : IComparable<T> { return MaxIndex(list, 0, list.Count); }
        public static Int2 MinMaxIndex<T>(this IList<T> list) where T : IComparable<T> { return MinMaxIndex(list, 0, list.Count); }

        public static int MinLastIndex<T>(this IList<T> list, int index, int length) where T : IComparable<T>
        {
            int j;
            for (int i = j = index + length - 1; --i >= index; )
                if (list[j].CompareTo(list[i]) > 0) j = i;
            return j;
        }
        public static int MaxLastIndex<T>(this IList<T> list, int index, int length) where T : IComparable<T>
        {
            int j;
            for (int i = j = index + length - 1; --i >= index; )
                if (list[j].CompareTo(list[i]) < 0) j = i;
            return j;
        }
        public static Int2 MinMaxLastIndex<T>(this IList<T> list, int index, int length) where T : IComparable<T>
        {
            int j0, j1;
            for (int i = j0 = j1 = index + length - 1; --i >= index; )
                if (list[j0].CompareTo(list[i]) > 0) j0 = i;
                else if (list[j1].CompareTo(list[i]) < 0) j1 = i;
            return new Int2(j0, j1);
        }
        public static int MinLastIndex<T>(this IList<T> list) where T : IComparable<T> { return MinLastIndex(list, 0, list.Count); }
        public static int MaxLastIndex<T>(this IList<T> list) where T : IComparable<T> { return MaxLastIndex(list, 0, list.Count); }
        public static Int2 MinMaxLastIndex<T>(this IList<T> list) where T : IComparable<T> { return MinMaxLastIndex(list, 0, list.Count); }

        public static int FindIndex<T>(this IList<T> list, Predicate<T> match)
        {
            for (int i = 0; i < list.Count; i++)
                if (match(list[i])) return i;
            return -1;
        }
        public static int FindIndexLast<T>(this IList<T> list, Predicate<T> match)
        {
            for (int i = list.Count; --i >= 0; )
                if (match(list[i])) return i;
            return -1;
        }
        public static int[] FindAllIndex<T>(this IList<T> list, Predicate<T> match)
        {
            var L = new List<int>();
            for (int i = 0; i < list.Count; i++)
                if (match(list[i])) L.Add(i);
            return L.ToArray();
        }

        // IList Others
        public static void Let<T>(this IList<T> list, T item)
        {
            for (int i = list.Count; --i >= 0; )
                list[i] = item;
        }
        public static T[] Sub<T>(this IList<T> list, int start)
        {
            T[] data = new T[list.Count - start];
            for (int i = 0; i < data.Length; i++)
                data[i] = list[i + start];
            return data;
        }
        public static T[] Sub<T>(this IList<T> list, int start, int count)
        {
            T[] data = new T[count];
            for (int i = 0; i < count; i++)
                data[i] = list[i + start];
            return data;
        }
        public static bool AllBackward<T>(this IList<T> list, Predicate<T> match)
        {
            for (int i = list.Count; --i >= 0; )
                if (!match(list[i])) return false;
            return true;
        }
        public static bool AnyBackward<T>(this IList<T> list, Predicate<T> match)
        {
            for (int i = list.Count; --i >= 0; )
                if (match(list[i])) return false;
            return true;
        }

        public static int BinarySearch<T>(this IList<T> list, T item) where T : IComparable<T>
        {
            int i0 = 0;
            int i1 = list.Count;
            while (true)
            {
                if (i0 == i1) return ~i0;
                int ii = (i0 + i1) / 2;
                int c = item.CompareTo(list[ii]);
                if (c == 0) return ii;
                if (c < 0) i1 = ii; else i0 = ii + 1;
            }
        }
        public static int BinarySearch<T>(this IList<T> list, T item, IComparer<T> comparer)
        {
            int i0 = 0;
            int i1 = list.Count;
            while (true)
            {
                if (i0 == i1) return ~i0;
                int ii = (i0 + i1) / 2;
                int c = comparer.Compare(item, list[ii]);
                if (c == 0) return ii;
                if (c < 0) i1 = ii; else i0 = ii + 1;
            }
        }

        public static T First<T>(this IList<T> list)
        {
            return list[0];
        }
        public static T Last<T>(this IList<T> list)
        {
            return list[list.Count - 1];
        }
        public static T Pop<T>(this IList<T> list)
        {
            T item = list[list.Count - 1];
            list.RemoveAt(list.Count - 1);
            return item;
        }

        public static int LetDistinctSorted<T>(this IList<T> list)
        {
            var comparer = EqualityComparer<T>.Default;
            int i = 0;
            for (int j = 1; j < list.Count; j++)
            {
                if (!comparer.Equals(list[i], list[j]))
                    list[++i] = list[j];
            }
            i++;
            int overlap = list.Count - i;
            for (int j = list.Count; --j >= i; )
                list.RemoveAt(j);
            return overlap;
        }

        public static void Overwrite<TSource>(this IList<TSource> list, Func<TSource, TSource> selector)
        {
            for (int i = 0; i < list.Count; i++)
                list[i] = selector(list[i]);
        }
        public static void Overwrite<TSource>(this IList<TSource> list, Func<TSource, int, TSource> selector)
        {
            for (int i = 0; i < list.Count; i++)
                list[i] = selector(list[i], i);
        }

        public static T[] Concatenate<T>(this IList<T> list, IList<T> other)
        {
            var array = new T[list.Count + other.Count];
            list.CopyTo(array, 0);
            other.CopyTo(array, list.Count);
            return array;
        }
        public static T[] Concatenate<T>(this IList<T> list, T item)
        {
            var array = new T[list.Count + 1];
            list.CopyTo(array, 0);
            array[list.Count] = item;
            return array;
        }
        public static T[] Concatenate<T>(T item, IList<T> list)
        {
            var array = new T[list.Count + 1];
            array[0] = item;
            list.CopyTo(array, 1);
            return array;
        }
        #endregion

        #region ICollection, SortedSet, IDictionary, ISet
        // ICollection
        public static void AddNew<T>(this ICollection<T> source) where T : new()
        {
            source.Add(new T());
        }
        public static void AddNew<T>(this ICollection<T> source, int count) where T : new()
        {
            for (int i = count; --i >= 0; )
                source.Add(new T());
        }

        // SortedSet
        public static T First<T>(this SortedSet<T> set)
        {
            return set.Min;
        }
        public static T Last<T>(this SortedSet<T> set)
        {
            return set.Max;
        }
        public static T Pop<T>(this SortedSet<T> set)
        {
            var item = set.Max;
            set.Remove(item);
            return item;
        }

        // IDictionary
        public static KeyValuePair<TKey, TValue> Pop<TKey, TValue>(this IDictionary<TKey, TValue> dictionary)
        {
            var item = dictionary.Last();
            dictionary.Remove(item.Key);
            return item;
        }
        public static bool AddOrDiscard<TKey, TValue>(this IDictionary<TKey, TValue> dictionary, TKey key, TValue value)
        {
            if (dictionary.ContainsKey(key)) return false;
            dictionary.Add(key, value); return true;
        }
        public static bool AddOrOverwrite<TKey, TValue>(this IDictionary<TKey, TValue> dictionary, TKey key, TValue value)
        {
            if (dictionary.ContainsKey(key)) { dictionary[key] = value; return false; }
            dictionary.Add(key, value); return true;
        }
        public static TValue GetItemOrDefault<TKey, TValue>(this IDictionary<TKey, TValue> dictionary, TKey key)
        {
            return dictionary.Keys.Contains(key) ? dictionary[key] : default(TValue);
        }
        public static TValue TryGetValue<TKey, TValue>(this IDictionary<TKey, TValue> dictionary, TKey key, TValue defaultvalue)
        {
            TValue value;
            return dictionary.TryGetValue(key, out value) ? value : defaultvalue;
        }
        public static TValue TryGetValue<TKey, TValue>(this IDictionary<TKey, TValue> dictionary, TKey key)
        {
            TValue value;
            return dictionary.TryGetValue(key, out value) ? value : default(TValue);
        }

        // ISet
        public static bool AddOrDiscard<T>(this ISet<T> set, T item)
        {
            return set.Add(item);
        }
        public static bool AddOrOverwrite<T>(this ISet<T> set, T item)
        {
            if (set.Contains(item)) { set.Remove(item); set.Add(item); return false; }
            else { set.Add(item); return true; }
        }

        // Dictionary
        // SortedDictionary
        // SortedList
        #endregion

        #region Sorting
        public static void Sort<T>(this T[] array) where T : IComparable<T> { Array.Sort(array); }
        public static void Sort<T>(this T[] array, bool ad) where T : IComparable<T> { if (ad) Array.Sort(array); else Array.Sort(array, (x, y) => y.CompareTo(x)); }
        public static void Sort<T>(this T[] array, Comparison<T> compare) { Array.Sort(array, compare); }
        public static void Sort<T>(this T[] array, Comparison<T> compare, bool ad) { if (ad) Array.Sort(array, compare); else Array.Sort(array, (x, y) => compare(y, x)); }

        public static void Sort<T>(this List<T> list, bool ad) where T : IComparable<T> { if (ad) list.Sort(); else list.Sort((x, y) => y.CompareTo(x)); }
        public static void Sort<T>(this List<T> list, Comparison<T> compare, bool ad) { list.Sort(ad ? compare : (x, y) => compare(y, x)); }
        public static void Sort<T>(this List<T> list, Func<T, int> selector) { list.Sort((x, y) => selector(x).CompareTo(selector(y))); }
        public static void Sort<T>(this List<T> list, Func<T, double> selector) { list.Sort((x, y) => selector(x).CompareTo(selector(y))); }
        public static void Sort<T>(this List<T> list, Func<T, int> selector, bool ad) { if (ad) list.Sort((x, y) => selector(x).CompareTo(selector(y))); else list.Sort((x, y) => selector(y).CompareTo(selector(x))); }
        public static void Sort<T>(this List<T> list, Func<T, double> selector, bool ad) { if (ad) list.Sort((x, y) => selector(x).CompareTo(selector(y))); else list.Sort((x, y) => selector(y).CompareTo(selector(x))); }

        public static int[] SortIndex<T>(this IList<T> list, Comparison<T> comparison)
        {
            int[] index = new int[list.Count];
            for (int i = index.Length; --i >= 0; )
                index[i] = i;
            Array.Sort(index, (x, y) => comparison(list[x], list[y]));
            return index;
        }
        public static int[] SortIndex<T>(this IList<T> list, Comparison<T> comparison, bool ad) { return ad ? SortIndex(list, comparison) : SortIndex(list, (x, y) => comparison(y, x)); }
        public static int[] SortIndex<T>(this IList<T> list) where T : IComparable<T> { return SortIndex(list, (x, y) => x.CompareTo(y)); }
        public static int[] SortIndex<T>(this IList<T> list, bool ad) where T : IComparable<T> { return ad ? SortIndex(list, (x, y) => x.CompareTo(y)) : SortIndex(list, (y, x) => x.CompareTo(y)); }
        public static int[] SortIndex<T>(this IList<T> list, Func<T, int> selector) { return SortIndex(list, (x, y) => selector(x).CompareTo(selector(y))); }
        public static int[] SortIndex<T>(this IList<T> list, Func<T, double> selector) { return SortIndex(list, (x, y) => selector(x).CompareTo(selector(y))); }
        public static int[] SortIndex<T>(this IList<T> list, Func<T, int> selector, bool ad) { return ad ? SortIndex(list, (x, y) => selector(x).CompareTo(selector(y))) : SortIndex(list, (x, y) => selector(y).CompareTo(selector(x))); }
        public static int[] SortIndex<T>(this IList<T> list, Func<T, double> selector, bool ad) { return ad ? SortIndex(list, (x, y) => selector(x).CompareTo(selector(y))) : SortIndex(list, (x, y) => selector(y).CompareTo(selector(x))); }

        public static int[] IndexToRank(IList<int> index)
        {
            int[] Rank = new int[index.Count];
            for (int i = index.Count; --i >= 0; )
                Rank[index[i]] = i;
            return Rank;
        }
        public static void SortAsIndex<T>(this IList<T> list, IList<int> index)
        {
            for (int i = index.Count; --i >= 0; )
            {
                if (index[i] < 0) continue;
                int k = i;
                T datum = list[k];
                while (true)
                {
                    int l = index[k]; index[k] = ~index[k];
                    if (l == i) break;
                    list[k] = list[l]; k = l;
                }
                list[k] = datum;
            }
            for (int i = index.Count; --i >= 0; )
                index[i] = ~index[i];
        }
        public static int Compare<T>(IList<T> left, IList<T> right) where T : IComparable<T>
        {
            if (left == null) if (right == null) return 0; else return -1; else if (right == null) return 1;
            if (left.Count != right.Count) return left.Count - right.Count;
            for (int i = 0; i < left.Count; i++)
            {
                int c = left[i].CompareTo(right[i]);
                if (c != 0) return c;
            }
            return 0;
        }
        public static int Compare<T>(T[] left, T[] right) where T : IComparable<T>
        {
            if (left == null) if (right == null) return 0; else return -1; else if (right == null) return 1;
            if (left.Length != right.Length) return left.Length - right.Length;
            for (int i = 0; i < left.Length; i++)
            {
                int c = left[i].CompareTo(right[i]);
                if (c != 0) return c;
            }
            return 0;
        }
        public static int Compare<T>(T[][] left, T[][] right) where T : IComparable<T>
        {
            if (left == null) if (right == null) return 0; else return -1; else if (right == null) return 1;
            if (left.Length != right.Length) return left.Length - right.Length;
            for (int i = 0; i < left.Length; i++)
            {
                int c = Compare(left[i], right[i]);
                if (c != 0) return c;
            }
            return 0;
        }
        public static int Compare<T>(T[][][] left, T[][][] right) where T : IComparable<T>
        {
            if (left == null) if (right == null) return 0; else return -1; else if (right == null) return 1;
            if (left.Length != right.Length) return left.Length - right.Length;
            for (int i = 0; i < left.Length; i++)
            {
                int c = Compare(left[i], right[i]);
                if (c != 0) return c;
            }
            return 0;
        }
        public static int Compare<T>(T[][][][] left, T[][][][] right) where T : IComparable<T>
        {
            if (left == null) if (right == null) return 0; else return -1; else if (right == null) return 1;
            if (left.Length != right.Length) return left.Length - right.Length;
            for (int i = 0; i < left.Length; i++)
            {
                int c = Compare(left[i], right[i]);
                if (c != 0) return c;
            }
            return 0;
        }
        #endregion

        #region TimeSpan
        public static TimeSpan Multiply(this TimeSpan left, int right) { return TimeSpan.FromTicks(left.Ticks * right); }
        public static TimeSpan Multiply(this TimeSpan left, double right) { return TimeSpan.FromTicks((long)(left.Ticks * right)); }
        public static TimeSpan Divide(this TimeSpan left, int right) { return TimeSpan.FromTicks(left.Ticks / right); }
        public static TimeSpan Divide(this TimeSpan left, double right) { return TimeSpan.FromTicks((long)(left.Ticks / right)); }
        public static TimeSpan Remainder(this TimeSpan left, int right) { return TimeSpan.FromTicks(left.Ticks % right); }
        public static TimeSpan Remainder(this TimeSpan left, double right) { return TimeSpan.FromTicks((long)(left.Ticks % right)); }
        public static TimeSpan Sum<TSource>(this IEnumerable<TSource> source, Func<TSource, TimeSpan> selector)
        {
            TimeSpan r = new TimeSpan();
            foreach (var element in source) checked { r += selector(element); }
            return r;
        }
        public static TimeSpan Average<TSource>(this IEnumerable<TSource> source, Func<TSource, TimeSpan> selector)
        {
            TimeSpan r = new TimeSpan();
            int count = 0;
            foreach (var element in source) { checked { r += selector(element); } count++; }
            return r.Divide(count);
        }
        #endregion

        #region string
        public static int TryParse(this string str, int defaultvalue)
        {
            int value;
            return int.TryParse(str, out value) ? value : defaultvalue;
        }
        public static double TryParse(this string str, double defaultvalue)
        {
            double value;
            return double.TryParse(str, out value) ? value : defaultvalue;
        }

        public static string Join(string separator, string terminator, IEnumerable<string> source)
        {
            var buffer = new StringBuilder();
            foreach (var element in source)
            {
                buffer.Append(element);
                buffer.Append(separator);
            }
            if (buffer.Length > 0)
                buffer.Remove(buffer.Length - separator.Length, separator.Length);
            buffer.Append(terminator);
            return buffer.ToString();
        }

        public static string[] Split(this string str, char separator) { return str.Split(new char[] { separator }); }
        public static string[] Split(this string str, string separator) { return str.Split(new string[] { separator }, StringSplitOptions.None); }

        public static string[] ToLines(this string str) { return str.Split(Gb.NewLineCodes, StringSplitOptions.None); }
        public static string[] ToItemsTab(this string line) { return line.Split('\t'); }
        public static string[] ToItemsCsv(this string line)
        {
            if (line == null) return null;
            var converted = new List<char>();
            bool quotation = false;
            bool escape = false;
            for (int i = 0; i < line.Length; i++)
            {
                char c = i < line.Length ? line[i] : ',';
                if (escape)
                {
                    if (c != '\"') { quotation = false; escape = false; }
                }
                if (!quotation)
                {
                    if (c == '\"') { quotation = true; continue; }
                    if (c == ',') c = '\t';
                }
                else
                {
                    if (c == '\"')
                    {
                        if (!escape) { escape = true; continue; }  // 1個目
                        else { escape = false; }                   // 2個目
                    }
                }
                converted.Add(c);
            }
            return new string(converted.ToArray()).ToItemsTab();
        }
        public static string ToLineTab(this IEnumerable<string> items) { return string.Join("\t", items); }
        public static string ToLineSpace(this IEnumerable<string> items) { return string.Join(" ", items); }

        public static string Replace(this string str, string[,] replacelist)
        {
            for (int i = 0; i < replacelist.GetLength(0); i++)
                str = str.Replace(replacelist[i, 0], replacelist[i, 1]);
            return str;
        }

        public static void Write(this StringBuilder sb, string value) { sb.Append(value); }
        public static void Write(this StringBuilder sb, bool value) { sb.Append(value.ToString()); }
        public static void Write(this StringBuilder sb, char value) { sb.Append(value.ToString()); }
        public static void Write(this StringBuilder sb, char[] buffer) { sb.Append(buffer, 0, buffer.Length); }
        public static void Write(this StringBuilder sb, char[] buffer, int index, int count) { sb.Append(buffer, index, count); }
        public static void Write(this StringBuilder sb, decimal value) { sb.Append(value.ToString()); }
        public static void Write(this StringBuilder sb, double value) { sb.Append(value.ToString()); }
        public static void Write(this StringBuilder sb, float value) { sb.Append(value.ToString()); }
        public static void Write(this StringBuilder sb, int value) { sb.Append(value.ToString()); }
        public static void Write(this StringBuilder sb, long value) { sb.Append(value.ToString()); }
        public static void Write(this StringBuilder sb, object value) { sb.Append(value.ToString()); }
        public static void Write(this StringBuilder sb, string format, object arg0) { sb.Append(string.Format(format, arg0)); }
        public static void Write(this StringBuilder sb, string format, object arg0, object arg1) { sb.Append(string.Format(format, arg0, arg1)); }
        public static void Write(this StringBuilder sb, string format, object arg0, object arg1, object arg2) { sb.Append(string.Format(format, arg0, arg1, arg2)); }
        public static void Write(this StringBuilder sb, string format, params object[] arg) { sb.Append(string.Format(format, arg)); }
        public static void Write(this StringBuilder sb, uint value) { sb.Append(value.ToString()); }
        public static void Write(this StringBuilder sb, ulong value) { sb.Append(value.ToString()); }

        public static void WriteLine(this StringBuilder sb) { sb.AppendLine(); }
        public static void WriteLine(this StringBuilder sb, string value) { sb.AppendLine(value); }
        public static void WriteLine(this StringBuilder sb, bool value) { sb.AppendLine(value.ToString()); }
        public static void WriteLine(this StringBuilder sb, char value) { sb.AppendLine(value.ToString()); }
        public static void WriteLine(this StringBuilder sb, char[] buffer) { sb.Append(buffer); sb.AppendLine(); }
        public static void WriteLine(this StringBuilder sb, char[] buffer, int index, int count) { sb.Append(buffer, index, count); sb.AppendLine(); }
        public static void WriteLine(this StringBuilder sb, decimal value) { sb.AppendLine(value.ToString()); }
        public static void WriteLine(this StringBuilder sb, double value) { sb.AppendLine(value.ToString()); }
        public static void WriteLine(this StringBuilder sb, float value) { sb.AppendLine(value.ToString()); }
        public static void WriteLine(this StringBuilder sb, int value) { sb.AppendLine(value.ToString()); }
        public static void WriteLine(this StringBuilder sb, long value) { sb.AppendLine(value.ToString()); }
        public static void WriteLine(this StringBuilder sb, object value) { sb.AppendLine(value.ToString()); }
        public static void WriteLine(this StringBuilder sb, string format, object arg0) { sb.AppendLine(string.Format(format, arg0)); }
        public static void WriteLine(this StringBuilder sb, string format, object arg0, object arg1) { sb.AppendLine(string.Format(format, arg0, arg1)); }
        public static void WriteLine(this StringBuilder sb, string format, object arg0, object arg1, object arg2) { sb.AppendLine(string.Format(format, arg0, arg1, arg2)); }
        public static void WriteLine(this StringBuilder sb, string format, params object[] arg) { sb.AppendLine(string.Format(format, arg)); }
        public static void WriteLine(this StringBuilder sb, uint value) { sb.AppendLine(value.ToString()); }
        public static void WriteLine(this StringBuilder sb, ulong value) { sb.AppendLine(value.ToString()); }
        #endregion

        #region other classes
        // Path
        public static string[] GetFiles(string path) { return GetFiles(path, SearchOption.TopDirectoryOnly); }
        public static string[] GetFiles(string path, SearchOption option)
        {
            string dirname = Path.GetDirectoryName(path);
            string filename = Path.GetFileName(path);
            if (filename.Length == 0) filename = "*";
            string[] paths = Directory.GetFiles(dirname, filename, option);
            return paths;
        }

        public static void Write(this Stream stream, byte[] array)
        {
            stream.Write(array, 0, array.Length);
        }
        public static void Write(this Stream stream, IEnumerable<byte> source)
        {
            var array = source.ToArray();
            stream.Write(array, 0, array.Length);
        }
        public static int Read(this FileStream filestream, byte[] array)
        {
            return filestream.Read(array, 0, array.Length);
        }

        // Average
        public static Double2 Average<TSource>(this IEnumerable<TSource> source, Func<TSource, Double2> selector)
        {
            Double2 r = default(Double2);
            int count = 0;
            foreach (var element in source) { checked { r += selector(element); } count++; }
            return r / count;
        }
        public static Double3 Average<TSource>(this IEnumerable<TSource> source, Func<TSource, Double3> selector)
        {
            Double3 r = default(Double3);
            int count = 0;
            foreach (var element in source) { checked { r += selector(element); } count++; }
            return r / count;
        }
        #endregion
    }
    #endregion

    #region Gb
    public static partial class Gb
    {
        public static readonly string NewLine = Environment.NewLine;
        public static readonly string[] NewLineCodes = { "\r\n", "\n", "\r" };
        public static readonly Encoding CodeUTF8 = Encoding.GetEncoding("utf-8");
        public static readonly Encoding CodeSJIS = Encoding.GetEncoding("shift_jis");
        public static readonly Encoding CodeEUC = Encoding.GetEncoding("euc-jp");
        public static readonly Encoding CodeJIS = Encoding.GetEncoding("iso-2022-jp");

        public static IEnumerable<string> EnumerateFileLines(string path, Encoding encoding)
        {
            using (var file = new StreamReader(path, encoding))
            {
                while (!file.EndOfStream)
                    yield return file.ReadLine();
            }
        }
        public static IEnumerable<string> EnumerateFileLines(string path)
        {
            return EnumerateFileLines(path, DetectEncoding(path));
        }
        public static Encoding DetectEncoding(string path)
        {
            var names = new string[] { "utf-8", "shift_jis", "euc-jp", "iso-2022-jp" };
            foreach (var name in names)
            {
                var encoding = Encoding.GetEncoding(name);
                bool sw = true;
                using (var reader = new StreamReader(path, encoding))
                {
                    while (!reader.EndOfStream)
                    {
                        string line = reader.ReadLine();
                        if (line.IndexOf('\ufffd') != -1) { sw = false; break; }
                    }
                    if (!sw) continue;
                    return encoding;
                }
            }
            //Console.Error.WriteLine("cannot detect proper encoding.");
            return null;
        }

        public static T FromString<T>(string str)
        {
            return (T)System.ComponentModel.TypeDescriptor.GetConverter(typeof(T)).ConvertFrom(str);
        }

        public static int[] ExtractNumbers(string str)
        {
            var list = new List<int>();
            string[] items = str.Split(',');
            foreach (string s in items)
            {
                if (s.IndexOf('-') < 0) { list.Add(int.Parse(s)); continue; }
                string[] t = s.Split('-');
                for (int i = int.Parse(t[0]), j = int.Parse(t[1]); i <= j; ++i)
                    list.Add(i);
            }
            list.Sort();
            return list.ToArray();
        }

    }

    #endregion

    #region StopWatch
    public class StopWatch
    {
        DateTime RegisteredTime;
        TimeSpan _Duration;
        bool _Running;
        public StopWatch() { Reset(); }

        public TimeSpan Duration { get { return _Running ? _Duration + (DateTime.Now - RegisteredTime) : _Duration; } }
        public bool Running { get { return _Running; } }

        public void Register() { RegisteredTime = DateTime.Now; }
        public void Reset() { _Duration = TimeSpan.Zero; Start(); }
        public void Start() { _Running = true; RegisteredTime = DateTime.Now; }
        public TimeSpan Stop() { _Running = false; return _Duration += DateTime.Now - RegisteredTime; }
        public override string ToString() { return Duration.ToString(); }
    }
    #endregion

    #region RandomMT
    public class RandomMT
    {
        const int N = 624, M = 397;
        const uint UPPER_MASK = 0x80000000u;
        const uint LOWER_MASK = 0x7fffffffu;
        const uint MATRIX_A = 0x9908b0dfu;

        uint BufferIndex;
        uint[] Buffer = new uint[N];

        public RandomMT() { Init((uint)DateTime.Now.ToBinary()); }
        public RandomMT(uint seed) { Init(seed); }
        public RandomMT(int seed) { Init((uint)seed); }

        public void Init(uint seed)
        {
            Buffer[0] = seed;
            for (int i = 1; i < N; i++)
                Buffer[i] = 1812433253u * (Buffer[i - 1] ^ (Buffer[i - 1] >> 30)) + (uint)i;
            BufferIndex = N;
        }
        public void Init(uint[] key)
        {
            if (key == null) ThrowException.ArgumentNullException("key");
            Init(19650218u);
            uint i = 1, j = 0;
            for (int k = Math.Max(N, key.Length); k > 0; k--)
            {
                Buffer[i] = (Buffer[i] ^ ((Buffer[i - 1] ^ (Buffer[i - 1] >> 30)) * 1664525u)) + key[j] + j; // non linear
                if (++i >= N) { i = 1; Buffer[0] = Buffer[N - 1]; }
                if (++j >= key.Length) j = 0;
            }
            for (int k = N - 1; k > 0; k--)
            {
                Buffer[i] = (Buffer[i] ^ ((Buffer[i - 1] ^ (Buffer[i - 1] >> 30)) * 1566083941u)) - i; // non linear
                if (++i >= N) { i = 1; Buffer[0] = Buffer[N - 1]; }
            }
            Buffer[0] = 0x80000000u; // MSB is 1; assuring non-zero initial array
        }

        // 32-bit uint [0, 0xffffffff]
        public uint UInt32()
        {
            if (BufferIndex >= N) { BufferIndex = 0; UInt32_(); } // generate N words at one time            
            uint y = Buffer[BufferIndex++];
            // tempering
            y ^= (y >> 11);
            y ^= (y << 7) & 0x9d2c5680u;
            y ^= (y << 15) & 0xefc60000u;
            y ^= (y >> 18);
            return y;
        }
        void UInt32_()
        {
            for (int kk = 0, k1 = 1, kM = M; kk < N; kk++, k1++, kM++)
            {
                if (k1 == N) k1 = 0;
                if (kM == N) kM = 0;
                uint y = (Buffer[kk] & UPPER_MASK) | (Buffer[k1] & LOWER_MASK);
                Buffer[kk] = Buffer[kM] ^ (y >> 1) ^ ((y & 1u) * MATRIX_A);
            }
        }
        // 31-bit int [-0x80000000, 0x7fffffff]
        public int Int32() { return (int)this.UInt32(); }
        // 31-bit int [0, 0x7fffffff]
        public int Int31() { return (int)(this.UInt32() >> 1); }

        public uint UInt(uint value) { return this.UInt32() % value; }
        public int Int(int value) { return (int)UInt((uint)value); }

        public uint UIntExact(uint value)
        {
            while (true)
            {
                uint y = this.UInt32();
                if (y < (1LU << 32) / value * value) return y % value;
            }
        }
        public int IntExact(int value) { return (int)UIntExact((uint)value); }

        public double Double() { return Double32CO(); }
        // 32-bit double [0, 1]
        public double Double32CC() { return this.UInt32() * (1.0 / ((1L << 32) - 1)); }
        // 32-bit double [0, 1)
        public double Double32CO() { return this.UInt32() * (1.0 / (1L << 32)); }
        // 32-bit double (0, 1]
        public double Double32OC() { return (this.UInt32() + 1.0) * (1.0 / (1L << 32)); }
        // 32-bit double (0, 1)
        public double Double32OO() { return (this.UInt32() + 0.5) * (1.0 / (1L << 32)); }
        // 53-bit double [0, 1]
        public double Double53CC() { return ((this.UInt32() >> 5) * (double)(1 << 26) + (this.UInt32() >> 6)) * (1.0 / ((1L << 53) - 1)); }
        // 53-bit double [0, 1)
        public double Double53CO() { return ((this.UInt32() >> 5) * (double)(1 << 26) + (this.UInt32() >> 6)) * (1.0 / (1L << 53)); }
        // 53-bit double (0, 1]
        public double Double53OC() { return ((this.UInt32() >> 5) * (double)(1 << 26) + (this.UInt32() >> 6) + 1.0) * (1.0 / (1L << 53)); }
        // 52-bit double (0, 1)
        public double Double52OO() { return ((this.UInt32() >> 6) * (double)(1 << 26) + (this.UInt32() >> 6) + 0.5) * (1.0 / (1L << 52)); }


        public double Gaussian() { return Gaussian32(); }

        double BufferGaussian32 = double.PositiveInfinity;
        public double Gaussian32()
        {
            if (BufferGaussian32 != double.PositiveInfinity) { double x = BufferGaussian32; BufferGaussian32 = double.PositiveInfinity; return x; }
            double v1, v2, r;
            do
            {
                v1 = this.Double32OO() - 0.5;
                v2 = this.Double32OO() - 0.5;
                r = v1 * v1 + v2 * v2;
            } while (r >= 0.25 || r == 0);
            double f = Math.Sqrt(-2 * Math.Log(r * 4) / r);
            BufferGaussian32 = v2 * f;
            return v1 * f;
        }
        double BufferGaussian52 = double.PositiveInfinity;
        public double Gaussian52()
        {
            if (BufferGaussian52 != double.PositiveInfinity) { double x = BufferGaussian52; BufferGaussian52 = double.PositiveInfinity; return x; }
            double v1, v2, r;
            do
            {
                v1 = this.Double52OO() - 0.5;
                v2 = this.Double52OO() - 0.5;
                r = v1 * v1 + v2 * v2;
            } while (r >= 0.25 || r == 0);
            double f = Math.Sqrt(-2 * Math.Log(r * 4) / r);
            BufferGaussian52 = v2 * f;
            return v1 * f;
        }
    }
    #endregion

    #region Sorting classes
    public class ComparerArray<T> : IComparer<T[]>
        where T : IComparable<T>
    {
        public ComparerArray() { }
        #region IComparer<T[]> メンバー
        public int Compare(T[] x, T[] y) { return Ex.Compare(x, y); }
        #endregion
    }
    public class ComparerArray2<T> : IComparer<T[][]>
        where T : IComparable<T>
    {
        public ComparerArray2() { }
        #region IComparer<T[][]> メンバー
        public int Compare(T[][] x, T[][] y) { return Ex.Compare(x, y); }
        #endregion
    }
    public class ComparerArray3<T> : IComparer<T[][][]>
        where T : IComparable<T>
    {
        public ComparerArray3() { }
        #region IComparer<T[][][]> メンバー
        public int Compare(T[][][] x, T[][][] y) { return Ex.Compare(x, y); }
        #endregion
    }
    public class ComparerArray4<T> : IComparer<T[][][][]>
        where T : IComparable<T>
    {
        public ComparerArray4() { }
        #region IComparer<T[][][]> メンバー
        public int Compare(T[][][][] x, T[][][][] y) { return Ex.Compare(x, y); }
        #endregion
    }
    public class ComparerIList<T> : IComparer<IList<T>>
        where T : IComparable<T>
    {
        public ComparerIList() { }
        #region IComparer<IList<T>> メンバー
        public int Compare(IList<T> x, IList<T> y) { return Ex.Compare(x, y); }
        #endregion
    }
    public class ComparerComparison<T> : IComparer<T>
    {
        Comparison<T> comparer;
        public ComparerComparison(Comparison<T> comparison) { this.comparer = comparison; }
        #region IComparer<T> メンバ
        public int Compare(T x, T y) { return comparer(x, y); }
        #endregion
    }
    public class ComparerReverse<T> : IComparer<T>
    {
        IComparer<T> comparer;
        public ComparerReverse() { this.comparer = Comparer<T>.Default; }
        public ComparerReverse(IComparer<T> comparer) { this.comparer = comparer; }
        public ComparerReverse(Comparison<T> comparer) { this.comparer = new ComparerComparison<T>(comparer); }
        #region IComparer<T> メンバ
        public int Compare(T x, T y) { return comparer.Compare(y, x); }
        #endregion
    }
    #endregion

    #region SortedList
    public class SortedList<T> : List<T>
    {
        protected readonly IComparer<T> _Comparer;

        public SortedList()
            : base()
        {
            this._Comparer = Comparer<T>.Default;
        }
        public SortedList(IComparer<T> comparer)
            : base()
        {
            this._Comparer = comparer;
        }
        public SortedList(int capacity)
            : base(capacity)
        {
            this._Comparer = Comparer<T>.Default;
        }
        public SortedList(IComparer<T> comparer, int capacity)
            : base(capacity)
        {
            this._Comparer = comparer;
        }

        public IComparer<T> Comparer
        {
            get { return _Comparer; }
        }

        public bool AddOrDiscard(T item)
        {
            int index = IndexOf(item);
            if (index >= 0) return false;
            else { Insert(~index, item); return true; }
        }
        public bool AddOrOverwrite(T item)
        {
            int index = IndexOf(item);
            if (index >= 0) { this[index] = item; return false; }
            else { Insert(~index, item); return true; }
        }
        public T Pop()
        {
            T item = this[Count - 1];
            RemoveAt(Count - 1);
            return item;
        }

        public new int IndexOf(T item) { return BinarySearch(item); }
        public new void Add(T item)
        {
            int index = BinarySearch(item);
            if (index >= 0) throw new ArgumentException();
            Insert(~index, item);
        }
        public new bool Contains(T item) { return BinarySearch(item) >= 0; }
    }
    #endregion

    #region ListedList
    public class ListedList<T> : IList<T>
    {
        protected const int FixedSize = 1024, FixedHalfSize = FixedSize / 2;
        protected List<List<T>> _Items = new List<List<T>>();
        protected int _Count = 0;

        public ListedList()
        {
        }
        public ListedList(IEnumerable<T> collections)
        {
            foreach (var item in collections) Add(item);
        }

        protected Int2 DecomposeIndex(int index)
        {
            int j = index;
            for (int i = 0; i < _Items.Count; i++)
            {
                if (j < _Items[i].Count) return new Int2(i, j);
                j -= _Items[i].Count;
            }
            ThrowException.ArgumentOutOfRangeException("index");
            return default(Int2);
        }
        protected int ComposeIndex(Int2 indexes)
        {
            int index = indexes.Y;
            for (int i = indexes.X; --i >= 0; )
                index += _Items[i].Count;
            return index;
        }
        protected Int2 _IndexOf(T item)
        {
            for (int i = 0; i < _Items.Count; ++i)
            {
                List<T> items = _Items[i];
                int j = _Items[i].IndexOf(item);
                if (j >= 0) return new Int2(i, j);
            }
            return ~new Int2(0, 0);
        }
        protected Int2 _BinarySearch(T item) { return _BinarySearch(item, Comparer<T>.Default); }
        protected Int2 _BinarySearch(T item, IComparer<T> comparer)
        {
            int i0 = 0, i1 = _Items.Count;
            while (i0 < i1)
            {
                int i = (i0 + i1) / 2;
                int c = comparer.Compare(item, _Items[i][0]);
                if (c == 0) return new Int2(i, 0);
                if (c < 0) i1 = i; else i0 = i + 1;
            }
            if (i0 == 0) return ~new Int2(0, 0);
            int ii = i0 - 1;
            List<T> items = _Items[ii];
            int j0 = 0, j1 = items.Count;
            while (j0 < j1)
            {
                int j = (j0 + j1) / 2;
                int c = comparer.Compare(item, items[j]);
                if (c == 0) return new Int2(ii, j);
                if (c < 0) j1 = j; else j0 = j + 1;
            }
            return ~new Int2(ii, j0);
        }

        protected void _Insert(Int2 indexes, T item)
        {
            if (indexes.X == _Items.Count) _Items.Add(new List<T>(FixedSize));
            if (_Items[indexes.X].Count == FixedSize)
            {
                _Items.Insert(indexes.X + 1, new List<T>(FixedSize));
                _Items[indexes.X + 1].AddRange(_Items[indexes.X].GetRange(FixedHalfSize, FixedHalfSize));
                _Items[indexes.X].RemoveRange(FixedHalfSize, FixedHalfSize);
                if (indexes.Y >= FixedHalfSize) indexes -= new Int2(-1, FixedHalfSize);
            }
            _Items[indexes.X].Insert(indexes.Y, item);
            _Count++;
        }
        protected void _RemoveAt(Int2 indexes)
        {
            _Count--;
            _Items[indexes.X].RemoveAt(indexes.Y);
            if (_Items[indexes.X].Count == 0) _Items.RemoveAt(indexes.X);
        }

        public int RemoveAll(Predicate<T> match)
        {
            int count = 0;
            _Items.RemoveAll(items =>
            {
                count += items.RemoveAll(match);
                return items.Count == 0;
            });
            _Count -= count;
            return count;
        }
        public void CopyTo(T[] array)
        {
            int index = 0;
            foreach (List<T> items in _Items) { items.CopyTo(array, index); index += items.Count; }
        }
        public IEnumerable<T> Reverse()
        {
            for (int i = _Items.Count; --i >= 0; )
            {
                List<T> items = _Items[i];
                for (int j = items.Count; --j >= 0; )
                    yield return items[j];
            }
        }
        public override string ToString()
        {
            return string.Format("Count = {0}", _Count);
        }

        public T Pop()
        {
            _Count--;
            List<T> items = _Items[_Items.Count - 1];
            T item = items[items.Count - 1];
            items.RemoveAt(items.Count - 1);
            if (items.Count == 0) _Items.RemoveAt(_Items.Count - 1);
            return item;
        }

        public void Sort()
        {
            var newlist = new ListedList<T>();
            foreach (var list in _Items) { list.Sort(); list.Reverse(); }
            var head = _Items.Select(list => list.Count - 1).ToList();
            var index = Enumerable.Range(0, _Items.Count).ToList();
            var comparer = new ComparerComparison<int>((y, x) => Comparer<T>.Default.Compare(_Items[x][head[x]], _Items[y][head[y]]));
            index.Sort(comparer);
            while (_Items.Count > 0)
            {
                int most = index.Pop();
                newlist.Add(_Items[most][head[most]]);
                if (--head[most] >= 0)
                {
                    int insert = ~index.BinarySearch(most, comparer);
                    if (insert < 0) insert = ~insert;
                    index.Insert(insert, most);
                }
                else
                {
                    for (int i = index.Count; --i >= 0; )
                        if (index[i] > most) index[i]--;
                    head.RemoveAt(most);
                    _Items.RemoveAt(most);
                }
            }
            _Items = newlist._Items;
        }

        public int LetDistinct()
        {
            var comparer = EqualityComparer<T>.Default;
            int overlap = 0;
            for (int i = _Items.Count; --i >= 0; )
            {
                var list = _Items[i];
                for (int j = list.Count; --j > 0; )
                {
                    if (comparer.Equals(list[j], list[j - 1])) { list.RemoveAt(j); overlap++; }
                }
                if (i > 0)
                {
                    if (comparer.Equals(list[0], _Items[i - 1].Last())) { list.RemoveAt(0); overlap++; }
                }
                if (list.Count == 0) _Items.RemoveAt(i);
            }
            _Count -= overlap;
            return overlap;
        }

        #region IList<T> メンバ
        public int IndexOf(T item)
        {
            Int2 indices = _IndexOf(item);
            if (indices.X >= 0) return ComposeIndex(indices);
            else return ~ComposeIndex(~indices);
        }
        public void Insert(int index, T item) { _Insert(DecomposeIndex(index), item); }
        public void RemoveAt(int index) { _RemoveAt(DecomposeIndex(index)); }
        public T this[int index]
        {
            get { Int2 indices = DecomposeIndex(index); return _Items[indices.X][indices.Y]; }
            set { Int2 indices = DecomposeIndex(index); _Items[indices.X][indices.Y] = value; }
        }
        #endregion
        #region ICollection<T> メンバ
        public void Add(T item)
        {
            if (_Items.Count == 0 || _Items[_Items.Count - 1].Count == FixedSize) _Items.Add(new List<T>(FixedSize));
            _Items[_Items.Count - 1].Add(item);
            _Count++;
        }
        public void Clear()
        {
            _Count = 0;
            _Items.Clear();
        }
        public bool Contains(T item)
        {
            return _IndexOf(item).X >= 0;
        }
        public void CopyTo(T[] array, int arrayIndex)
        {
            int index = arrayIndex;
            foreach (var items in _Items)
            {
                items.CopyTo(array, index);
                index += items.Count;
            }
        }
        public int Count
        {
            get { return _Count; }
        }
        public bool IsReadOnly
        {
            get { return false; }
        }
        public bool Remove(T item)
        {
            Int2 indices = _IndexOf(item);
            if (indices.X >= 0) { _RemoveAt(indices); return true; }
            else { return false; }
        }
        #endregion
        #region IEnumerable<T> メンバ
        public IEnumerator<T> GetEnumerator()
        {
            for (int i = 0; i < _Items.Count; ++i)
            {
                List<T> items = _Items[i];
                for (int j = 0; j < items.Count; ++j)
                    yield return items[j];
            }
        }
        #endregion
        #region IEnumerable メンバ
        IEnumerator IEnumerable.GetEnumerator() { return GetEnumerator(); }
        #endregion
    }
    #endregion

    #region SortedListedList
    public class SortedListedList<T> : ListedList<T>
    {
        protected readonly IComparer<T> _Comparer;

        public SortedListedList()
            : base()
        {
            this._Comparer = Comparer<T>.Default;
        }
        public SortedListedList(IEnumerable<T> collections)
        {
            foreach (var item in collections) Add(item);
        }
        public SortedListedList(IComparer<T> comparer)
            : base()
        {
            this._Comparer = comparer;
        }

        public IComparer<T> Comparer
        {
            get { return _Comparer; }
        }

        protected new Int2 _BinarySearch(T item) { return _BinarySearch(item, _Comparer); }

        public int BinarySearch(T item) { return ComposeIndex(_BinarySearch(item)); }
        public new bool Add(T item)
        {
            Int2 indices = _BinarySearch(item);
            if (indices.X >= 0) { return false; }
            else { _Insert(~indices, item); return true; }
        }
        public new void Insert(int index, T item) { throw new NotImplementedException(); }

        public bool AddOrDiscard(T item)
        {
            Int2 indices = _BinarySearch(item);
            if (indices.X >= 0) return false;
            else { _Insert(~indices, item); return true; }
        }
        public bool AddOrOverwrite(T item)
        {
            Int2 indices = _BinarySearch(item);
            if (indices.X >= 0) { _Items[indices.X][indices.Y] = item; return false; }
            else { _Insert(~indices, item); return true; }
        }
        public T FindOrDefault(T item)
        {
            Int2 index = _BinarySearch(item);
            return index.X >= 0 ? _Items[index.X][index.Y] : default(T);
        }

        public new bool Contains(T item)
        {
            return _BinarySearch(item).X >= 0;
        }
    }
    #endregion

    #region WaveFile
    public static class WaveFile
    {
        static readonly char[] chunkRiff = { 'R', 'I', 'F', 'F' };
        static readonly char[] chunkType = { 'W', 'A', 'V', 'E' };
        static readonly char[] chunkFrmt = { 'f', 'm', 't', ' ' };
        static readonly char[] chunkData = { 'd', 'a', 't', 'a' };

        public static void Save(string filename, double[][] data, int samplesPerSecond, int bitsPerSample = 16)
        {
            try
            {
                using (var writer = new BinaryWriter(File.Open(filename, FileMode.Create)))
                {
                    int channels = data.Length;
                    int dataLength = data[0].Length;

                    int padding = 1;                           // File padding
                    int bytesPerSample = (bitsPerSample / 8) * channels;  // Bytes per sample.
                    int averageBytesPerSecond = bytesPerSample * samplesPerSecond;
                    int chunkDataLength = (bitsPerSample / 8) * dataLength * channels;
                    int chunkFrmtLength = 16;                    // Format chunk length.
                    int chunkRiffLength = chunkDataLength + 36;  // File length, minus first 8 bytes of RIFF description.

                    writer.Write(chunkRiff);                     //4 bytes
                    writer.Write(chunkRiffLength);               //4 bytes
                    {
                        writer.Write(chunkType);                     //4 bytes

                        writer.Write(chunkFrmt);                     //4 bytes
                        writer.Write(chunkFrmtLength);               //4 bytes
                        {
                            writer.Write((short)padding);                //2 bytes
                            writer.Write((short)channels);               //2 bytes
                            writer.Write(samplesPerSecond);              //4 bytes
                            writer.Write(averageBytesPerSecond);         //4 bytes
                            writer.Write((short)bytesPerSample);         //2 bytes
                            writer.Write((short)bitsPerSample);          //2 bytes
                        }
                        writer.Write(chunkData);                     //4 bytes
                        writer.Write(chunkDataLength);               //4 bytes
                        switch (bitsPerSample)
                        {
                            case 8: for (int i = 0; i < dataLength; ++i) for (int j = 0; j < channels; ++j) writer.Write((Byte)Mt.MinMax(data[j][i] * 0x80 + 0x80, 0, 0xff)); break;
                            case 16: for (int i = 0; i < dataLength; ++i) for (int j = 0; j < channels; ++j) writer.Write((Int16)Mt.MinMax(data[j][i] * 0x8000, -0x8000, 0x7fff)); break;
                            case 32: for (int i = 0; i < dataLength; ++i) for (int j = 0; j < channels; ++j) writer.Write((Int32)Mt.MinMax(data[j][i] * 0x80000000L, -0x80000000L, 0x7fffffffL)); break;
                            default: throw new Exception("WaveFormat.Save: unknown BitsPerSample");
                        }
                    }
                }
            }
            catch (Exception e)
            {
                Console.Error.WriteLine(e.Message);
            }
        }

        public static Tuple<double[][], int> Load(string filename)
        {
            try
            {
                using (var reader = new BinaryReader(File.Open(filename, FileMode.Open)))
                {
                    if (Ex.Compare(reader.ReadChars(4), chunkRiff) != 0) throw new Exception("WaveFormat.Load: not RIFF chunk");
                    reader.ReadInt32();  //chunkRiffLength
                    if (Ex.Compare(reader.ReadChars(4), chunkType) != 0) throw new Exception("WaveFormat.Load: not WAVE chunk");

                    int bitsPerSample = 0;
                    int samplesPerSecond = 0;
                    int channels = 0;
                    while (true)
                    {
                        char[] chunkName = reader.ReadChars(4);
                        int chunkLength = reader.ReadInt32();
                        if (Ex.Compare(chunkName, chunkFrmt) == 0)
                        {
                            reader.ReadInt16();  //shPad                 //2 bytes
                            channels = reader.ReadInt16();               //2 bytes
                            samplesPerSecond = reader.ReadInt32();       //4 bytes
                            reader.ReadInt32();  //averageBytesPerSecond //4 bytes
                            reader.ReadInt16();  //shBytesPerSample      //2 bytes
                            bitsPerSample = reader.ReadInt16();          //2 bytes
                            if (chunkLength > 16) reader.ReadBytes(chunkLength - 16);  //unknown data
                        }
                        else if (Ex.Compare(chunkName, chunkData) == 0)
                        {
                            int dataLength = chunkLength / channels / (bitsPerSample / 8);
                            var data = New.Array(channels, i => new double[dataLength]);
                            switch (bitsPerSample)
                            {
                                case 8: for (int i = 0; i < dataLength; ++i) for (int j = 0; j < channels; ++j) data[j][i] = (reader.ReadByte() - 0x80) / (double)0x80; break;
                                case 16: for (int i = 0; i < dataLength; ++i) for (int j = 0; j < channels; ++j) data[j][i] = reader.ReadInt16() / (double)0x8000; break;
                                case 32: for (int i = 0; i < dataLength; ++i) for (int j = 0; j < channels; ++j) data[j][i] = reader.ReadInt32() / (double)0x80000000L; break;
                                default: throw new Exception("WaveFormat.Load: unknown BitsPerSample");
                            }
                            return new Tuple<double[][], int>(data, samplesPerSecond);
                        }
                        else
                        {
                            reader.ReadBytes(chunkLength);  //unknown chunk
                        }
                    }
                }
            }
            catch (Exception e)
            {
                Console.Error.WriteLine(e.Message);
            }
            return null;
        }
    }
    #endregion

    public static partial class Mt
    {
        #region miscellaneous
        public const double DoubleEpsilon = 2.2250738585072013830902327173324e-308;  // 2^-1022
        public const double DoubleEps = 2.2204460492503130808472633361816e-16;  // 2^-52
        public const double DoubleFpMin = 1.0020841800044863889980540256751e-292;  // 2^-1022 / 2^-52 = 2^-970

        public static bool IsTooSmall(double x, double y) { return x + y == y; }

        public static void Swap<T>(ref T val1, ref T val2) { T z = val1; val1 = val2; val2 = z; }

        public static T MinMax<T>(T value, T valmin, T valmax) where T : IComparable<T> { return value.CompareTo(valmin) < 0 ? valmin : (value.CompareTo(valmax) > 0 ? valmax : value); }
        public static void LetMin<T>(ref T value, T valmin) where T : IComparable<T> { if (value.CompareTo(valmin) > 0) value = valmin; }
        public static void LetMax<T>(ref T value, T valmax) where T : IComparable<T> { if (value.CompareTo(valmax) < 0) value = valmax; }
        public static void LetMinMax<T>(ref T value, T valmin, T valmax) where T : IComparable<T> { if (value.CompareTo(valmin) < 0) value = valmin; else if (value.CompareTo(valmax) > 0) value = valmax; }
        public static bool IfLetMin<T>(ref T value, T valmin) where T : IComparable<T> { if (value.CompareTo(valmin) > 0) { value = valmin; return true; } return false; }
        public static bool IfLetMax<T>(ref T value, T valmax) where T : IComparable<T> { if (value.CompareTo(valmax) < 0) { value = valmax; return true; } return false; }

        public static Int2 DivRem(int left, int right) { int rem, div = Math.DivRem(left, right, out rem); return new Int2(rem, div); }
        public static int AlignUp(int value, int size) { return (value + (size - 1)) - (value + (size - 1)) % size; }
        public static int AlignDown(int value, int size) { return value - value % size; }
        public static int AlignUpX(int value, int size) { return (value + (size - 1)) & (~(size - 1)); }
        public static int AlignDownX(int value, int size) { return value & (~(size - 1)); }
        public static int RoundOff(double value) { return (int)Math.Floor(value + 0.5); }
        public static double RoundOff(double value, int digit) { double k = Math.Pow(10, digit); return Math.Floor(value * k + 0.5) / k; }
        public static Int2 RoundOff(Double2 value) { return new Int2(Mt.RoundOff(value.X), Mt.RoundOff(value.Y)); }
        public static Int3 RoundOff(Double3 value) { return new Int3(Mt.RoundOff(value.X), Mt.RoundOff(value.Y), Mt.RoundOff(value.Z)); }

        public static int MinMaxC(int value, int valmin, int valmax) { return value < valmin ? valmin : value >= valmax ? valmax - 1 : value; }
        public static void LetMinMaxC(ref int value, int valmin, int valmax) { if (value < valmin) value = valmin; else if (value >= valmax) value = valmax - 1; }
        public static int Div_(int left, int right) { return (left - ((left % right) + right) % right) / right; }  //right>0:下方向商   right<0:上方向商
        public static int Mod_(int left, int right) { return ((left % right) + right) % right; }        //right>0:下方向剰余 right<0:上方向剰余        
        //  5/ 3= 1	Div_( 5, 3)= 1
        // -5/ 3=-1	Div_(-5, 3)=-2
        //  5/-3=-1	Div_( 5,-3)=-2
        // -5/-3= 1	Div_(-5,-3)= 1
        //  5% 3= 2	Mod_( 5, 3)= 2
        // -5% 3=-2	Mod_(-5, 3)= 1
        //  5%-3= 2	Mod_( 5,-3)=-1
        // -5%-3=-2	Mod_(-5,-3)=-2

        public static int Sq(int value) { return value * value; }
        public static double Sq(double value) { return value * value; }
        public static int Cube(int value) { return value * value * value; }
        public static double Cube(double value) { return value * value * value; }
        public static int LetAbs(ref int value) { if (value < 0) value *= -1; return value; }
        public static double LetAbs(ref double value) { if (value < 0) value *= -1; return value; }
        public static void LetOrder(ref int val1, ref int val2) { if (val2 < val1) Swap(ref val1, ref val2); }
        public static void LetOrder(ref double val1, ref double val2) { if (val2 < val1) Swap(ref val1, ref val2); }

        public static int DivCeil(int left, int right) { return (left + right - 1) / right; }
        public static long DivCeil(long left, int right) { return (left + right - 1) / right; }
        public static int DivRound(int left, int right) { return (left + (left >= 0 ? right / 2 : -right / 2)) / right; }
        public static long DivRound(long left, int right) { return (left + (left >= 0 ? right / 2 : -right / 2)) / right; }

        public static int Min(int val1, int val2, int val3) { return Math.Min(Math.Min(val1, val2), val3); }
        public static int Max(int val1, int val2, int val3) { return Math.Max(Math.Max(val1, val2), val3); }
        public static int Min(int val1, int val2, int val3, int val4) { return Math.Min(Math.Min(val1, val2), Math.Min(val3, val4)); }
        public static int Max(int val1, int val2, int val3, int val4) { return Math.Max(Math.Max(val1, val2), Math.Max(val3, val4)); }
        public static double Min(double val1, double val2, double val3) { return Math.Min(Math.Min(val1, val2), val3); }
        public static double Max(double val1, double val2, double val3) { return Math.Max(Math.Max(val1, val2), val3); }
        public static double Min(double val1, double val2, double val3, double val4) { return Math.Min(Math.Min(val1, val2), Math.Min(val3, val4)); }
        public static double Max(double val1, double val2, double val3, double val4) { return Math.Max(Math.Max(val1, val2), Math.Max(val3, val4)); }

        public static double SlightlyInferior(double value) { return value - Math.Max(1e-296, Math.Abs(value) * 1e-12); }
        public static double SlightlySuperior(double value) { return value + Math.Max(1e-296, Math.Abs(value) * 1e-12); }
        public static int Log2Int(int value)
        {
            if (value <= 0) ThrowException.ArgumentOutOfRangeException("value");
            return Log2Int((uint)value);
        }
        public static int Log2Int(uint value)
        {
            if (value == 0) ThrowException.ArgumentOutOfRangeException("value");
            int i = 16;
            for (int j = 8; j > 0; j >>= 1)
                if (value < (1u << i)) i -= j; else i += j;
            return value < (1u << i) ? i - 1 : i;
        }
        #endregion

        //---------------------------------------------------------------------------
        #region Median functions
        public static double QuantileSorted(IList<double> data, double quantile)
        {
            int N = data.Count;
            if (N <= 0) return 0;
            if (N == 1) return data[0];
            double R = (N - 1) * MinMax(quantile, 0, 1);
            int M = (int)R;
            R -= M;
            return R == 0 ? data[M] : data[M] * (1 - R) + data[M + 1] * R;
        }

        public static double Quantile(IList<double> data, double quantile)
        {
            int N = data.Count;
            if (N <= 0) ThrowException.NoElements();
            if (N == 1) return data[0];

            double R = (N - 1) * MinMax(quantile, 0, 1);
            int M = (int)R, N0 = M + 1, N1 = N - N0;
            R -= M;
            bool I0, I1, F0, F1;
            I0 = I1 = F0 = F1 = false;
            double D0, D1;
            D0 = D1 = data[0];
            for (int j = N; --j >= 0; )
            {
                double D = data[j];
                if (I0 && (!(D0 < D)) || I1 && (!(D < D1))) continue;

                int LT = 0, GT = 0;
                for (int i = N - 1; ; )
                {
                    double d = data[i];
                    if (D < d) { if (++GT > N1) { I0 = true; D0 = D; break; } }
                    else
                        if (d < D) { if (++LT > N0) { I1 = true; D1 = D; break; } }
                    if (--i >= 0) continue;
                    if (!F0 && LT < N0) { I0 = F0 = true; D0 = D; if (F1) j = 0; }
                    if (!F1 && GT < N1) { I1 = F1 = true; D1 = D; if (F0) j = 0; }
                    break;
                }
            }
            return D0 * (1 - R) + D1 * R;
        }
        public static double Median(IList<double> data)
        {
            return Quantile(data, 0.5);
        }
        #endregion

        //---------------------------------------------------------------------------
        #region Triangular functions
        public static double Atanh(double value)
        {
            return (Math.Abs(value) >= 1 ? (value > 0 ? double.PositiveInfinity : double.NegativeInfinity) : Math.Log((1 + value) / (1 - value)) / 2);
        }
        public static double Atanh_(double value)
        {
            return (Math.Abs(value) >= 1 ? (value > 0 ? 18.7149738751185 : -18.7149738751185) : Math.Log((1 + value) / (1 - value)) / 2);
        }
        public static double Tanh(double value)
        {
            if (value == double.PositiveInfinity) return 1;
            if (value == double.NegativeInfinity) return -1;
            double y = Math.Exp(value * 2);
            return (y - 1) / (y + 1);
        }
        public static double Tanh_(double value)
        {
            if (Math.Abs(value) >= 18.7149738751185) return value > 0 ? 1 - 1e-16 : -1 + 1e-16;
            double y = Math.Exp(value + value);
            return (y - 1) / (y + 1);
        }
        public static double Acos_(double value)
        {
            return value <= -1 ? Math.PI : value >= 1 ? 0.0 : Math.Acos(value);
        }
        public static double Asin_(double value)
        {
            return value <= -1 ? -Math.PI / 2 : value >= 1 ? Math.PI / 2 : Math.Asin(value);
        }
        public static double Atan1(double y, double x)
        {
            double theta = Math.Atan2(y, x);
            if (theta <= -Math.PI / 2) return theta + Math.PI;
            if (theta > Math.PI / 2) return theta - Math.PI;
            return theta;
        }
        #endregion

        //---------------------------------------------------------------------------
        #region Special functions

        // ガンマ関数
        readonly static double[] LogGammaCoefficients = { 57.1562356658629235, -59.5979603554754912, 14.1360979747417471, -0.491913816097620199, .339946499848118887e-4, .465236289270485756e-4, -.983744753048795646e-4, .158088703224912494e-3, -.210264441724104883e-3, .217439618115212643e-3, -.164318106536763890e-3, .844182239838527433e-4, -.261908384015814087e-4, .368991826595316234e-5 };
        public static double LogGamma(double value)
        {
            if (value <= 0) return double.NaN;
            double ser = 0.999999999999997092;
            for (int j = 0; j < LogGammaCoefficients.Length; j++) ser += LogGammaCoefficients[j] / (j + 1 + value);
            double tmp = value + 5.24218750000000000;
            return (value + 0.5) * Math.Log(tmp) - tmp + Math.Log(2.5066282746310005 * ser / value);
        }
        public static double Gamma(double value) { return Math.Exp(LogGamma(value)); }

        // 階乗
        // 22! まではdoubleで正確に表現可能
        // 170! まではdoubleで近似的に表現可能
        static double[] FactorialBuffer;
        static void Factorial_()
        {
            FactorialBuffer = new double[171];
            double f = 1;
            FactorialBuffer[0] = f;
            for (int i = 1; i < FactorialBuffer.Length; i++)
            {
                f *= i;
                FactorialBuffer[i] = f;
            }
        }
        public static double Factorial(int value)
        {
            if (FactorialBuffer == null) Factorial_();
            if (value < 0 || value >= FactorialBuffer.Length) ThrowException.ArgumentOutOfRangeException("value");
            return FactorialBuffer[value];
        }

        static double[] LogFactorialBuffer;
        static void LogFactorial_()
        {
            LogFactorialBuffer = new double[2000];
            double f = 0;
            for (int i = 2; i < LogFactorialBuffer.Length; i++)
            {
                f += Math.Log(i);
                LogFactorialBuffer[i] = i <= 22 ? Math.Log(Factorial(i)) : f;
            }
        }
        public static double LogFactorial(int value)
        {
            if (value < 0) ThrowException.ArgumentOutOfRangeException("value");
            if (LogFactorialBuffer == null) LogFactorial_();
            if (value < LogFactorialBuffer.Length) return LogFactorialBuffer[value];
            return Mt.LogGamma(value + 1.0);
        }

        // Pochhammerと同じ
        public static double RisingFactorial(int value, int count)
        {
            if (count == 0) return 1;
            if (count == 1) return value;
            int s, v0, v1;
            if (value > 0)
            {
                s = 1;
                v0 = value + count - 1;
                v1 = value - 1;
            }
            else
            {
                s = 1 - 2 * (count & 1);
                v0 = -value;
                v1 = -value - count;
            }
            if (v0 < 0) return double.NaN;
            if (v1 < 0) return 0;
            return s * (v1 < 171 && v0 < 171 ? Factorial(v0) / Factorial(v1) : Math.Exp(LogFactorial(v0) - LogFactorial(v1)));
        }
        // FactorialPowerと同じ
        public static double FallingFactorial(int value, int count)
        {
            if (count == 0) return 1;
            if (count == 1) return value;
            int s, v0, v1;
            if (value >= count)
            {
                s = 1;
                v0 = value;
                v1 = value - count;
            }
            else
            {
                s = 1 - 2 * (count & 1);
                v0 = -value + count - 1;
                v1 = -value - 1;
            }
            if (v0 < 0) return double.NaN;
            if (v1 < 0) return 0;
            return s * (v1 < 171 && v0 < 171 ? Factorial(v0) / Factorial(v1) : Math.Exp(LogFactorial(v0) - LogFactorial(v1)));
        }

        public static double LogBinomial(int left, int right)
        {
            if (left < 0) ThrowException.ArgumentOutOfRangeException("left");
            int val2 = left - right;
            if (right < 0 || val2 < 0) return double.NegativeInfinity;
            if (right == 0 || val2 == 0) return 0.0;
            return Mt.LogFactorial(left) - (Mt.LogFactorial(right) + Mt.LogFactorial(val2));
        }
        public static double Binomial(int left, int right)
        {
            if (left < 0) ThrowException.ArgumentOutOfRangeException("left");
            int val2 = left - right;
            if (right < 0 || val2 < 0) return 0.0;
            if (right == 0 || val2 == 0) return 1.0;
            return Math.Round(Math.Exp(Mt.LogFactorial(left) - (Mt.LogFactorial(right) + Mt.LogFactorial(val2))), MidpointRounding.AwayFromZero);
        }

        public static double LogMultinomial(IEnumerable<int> source)
        {
            int total = 0;
            double sum = 0;
            foreach (int element in source)
            {
                if (element < 0) ThrowException.ArgumentOutOfRangeException("element");
                total += element;
                sum += Mt.LogFactorial(element);
            }
            return Mt.LogFactorial(total) - sum;
        }
        public static double Multinomial(IEnumerable<int> table)
        {
            return Math.Round(Math.Exp(LogMultinomial(table)), MidpointRounding.AwayFromZero);
        }

        // ベータ関数
        public static double LogBeta(double val1, double val2)
        {
            return Mt.LogGamma(val1) + Mt.LogGamma(val2) - Mt.LogGamma(val1 + val2);
        }
        public static double Beta(double val1, double val2) { return Math.Exp(LogBeta(val1, val2)); }


        public static double Logistic(double value)
        {
            return 1 / (1 + Math.Exp(-value));
        }

        // 不完全ガンマ関数
        readonly static double[] Gauleg18y = new double[] {
            0.0021695375159141994, 0.011413521097787704, 0.027972308950302116, 0.051727015600492421,
            0.082502225484340941, 0.12007019910960293, 0.16415283300752470, 0.21442376986779355,
            0.27051082840644336, 0.33199876341447887, 0.39843234186401943, 0.46931971407375483,
            0.54413605556657973, 0.62232745288031077, 0.70331500465597174, 0.78649910768313447,
            0.87126389619061517, 0.95698180152629142
        };
        readonly static double[] Gauleg18w = new double[] {
            0.0055657196642445571, 0.012915947284065419, 0.020181515297735382, 0.027298621498568734,
            0.034213810770299537, 0.040875750923643261, 0.047235083490265582, 0.053244713977759692,
            0.058860144245324798, 0.064039797355015485, 0.068745323835736408, 0.072941885005653087,
            0.076598410645870640, 0.079687828912071670, 0.082187266704339706, 0.084078218979661945,
            0.085346685739338721, 0.085983275670394821
        };
        static double gammpapprox(double gamma, double value, bool psig)
        {
            double a1 = gamma - 1, lna1 = Math.Log(a1), sqrta1 = Math.Sqrt(a1);
            double xu = (value > a1) ?
                Math.Max(a1 + 11.5 * sqrta1, value + 6.0 * sqrta1) :
                Math.Max(0.0, Math.Min(a1 - 7.5 * sqrta1, value - 5.0 * sqrta1));
            double sum = 0;
            for (int j = 0; j < Gauleg18y.Length; j++)
            {
                double t = value + (xu - value) * Gauleg18y[j];
                sum += Gauleg18w[j] * Math.Exp(a1 - t + a1 * (Math.Log(t) - lna1));
            }
            double ans = sum * (xu - value) * Math.Exp(a1 * (lna1 - 1.0) - LogGamma(gamma));
            return psig ? (ans > 0.0 ? 1.0 : 0.0) - ans : (ans >= 0.0 ? 0.0 : 1.0) + ans;
        }
        static double gser(double gamma, double value)
        {
            double sum = 1.0 / gamma;
            double del = sum;
            for (int i = 1; ; i++)
            {
                del *= value / (gamma + i);
                sum += del;
                if (Math.Abs(del) < Math.Abs(sum) * DoubleEps) break;
            }
            return sum * Math.Exp(gamma * Math.Log(value) - value - LogGamma(gamma));
        }
        static double gcf(double gamma, double value)
        {
            double c = double.PositiveInfinity;
            double d = 1 / (value - gamma + 1);
            double h = d;
            for (int n = 1; ; n++)
            {
                double an = n * (gamma - n);
                double bn = (value - gamma + 1) + n * 2;
                c = bn + an / c; if (Math.Abs(c) < DoubleFpMin) c = DoubleFpMin;
                d = bn + an * d; if (Math.Abs(d) < DoubleFpMin) d = DoubleFpMin;
                if (c == d) break;
                h *= c / d;
                d = 1.0 / d;
            }
            return h * Math.Exp(gamma * Math.Log(value) - value - LogGamma(gamma));
        }

        // = ∫0~value t^(gamma-1) e^-t dt / Γ(gamma)
        public static double IncompleteGammaLower(double gamma, double value)
        {
            if (gamma <= 0.0) return double.NaN;
            if (value == double.PositiveInfinity) return 1.0;
            if (value <= 0.0) return 0.0;
            if (gamma >= 100.0) return gammpapprox(gamma, value, true);
            return (value < gamma + 1.0) ? gser(gamma, value) : 1.0 - gcf(gamma, value);
        }
        // = ∫value~∞ t^(gamma-1) e^-t dt / Γ(gamma)
        public static double IncompleteGammaUpper(double gamma, double value)
        {
            if (gamma <= 0.0) return double.NaN;
            if (value == double.PositiveInfinity) return 0.0;
            if (value <= 0.0) return 1.0;
            if (gamma >= 100.0) return gammpapprox(gamma, value, false);
            return (value < gamma + 1.0) ? 1.0 - gser(gamma, value) : gcf(gamma, value);
        }
        static double InverseIncompleteGamma(double value, double gamma)
        {
            if (gamma <= 0.0) return double.NaN;
            if (value >= 1.0) return Math.Max(100.0, gamma + 100.0 * Math.Sqrt(gamma));
            if (value <= 0.0) return 0.0;

            double x, t, lna1 = 0, afac = 0, a1 = gamma - 1;
            const double EPS_ = 1.0e-8;
            double gln = Mt.LogGamma(gamma);
            if (gamma > 1.0)
            {
                lna1 = Math.Log(a1);
                afac = Math.Exp(a1 * (lna1 - 1.0) - gln);
                t = Math.Sqrt(-2.0 * Math.Log((value < 0.5) ? value : 1.0 - value));
                x = (2.30753 + t * 0.27061) / (1.0 + t * (0.99229 + t * 0.04481)) - t;
                if (value < 0.5) x = -x;
                x = Math.Max(1.0e-3, gamma * Math.Pow(1.0 - 1.0 / (9.0 * gamma) - x / (3.0 * Math.Sqrt(gamma)), 3));
            }
            else
            {
                t = 1.0 - gamma * (0.253 + gamma * 0.12);
                x = (value < t) ?
                    Math.Pow(value / t, 1.0 / gamma) :
                    1.0 - Math.Log(1.0 - (value - t) / (1.0 - t));
            }
            for (int j = 0; j < 12; j++)
            {
                if (x <= 0.0) return 0.0;
                double err = IncompleteGammaLower(gamma, x) - value;
                t = (gamma > 1.0) ?
                    afac * Math.Exp(-(x - a1) + a1 * (Math.Log(x) - lna1)) :
                    Math.Exp(-x + a1 * Math.Log(x) - gln);
                double u = err / t;
                x -= (t = u / (1.0 - 0.5 * Math.Min(1.0, u * (a1 / x - 1))));
                if (x <= 0.0) x = 0.5 * (x + t);
                if (Math.Abs(t) < EPS_ * x) break;
            }
            return x;
        }

        // ChiSquare分布
        public static double ChiSquareDistributionLower(double freedom, double value)
        {
            return Mt.IncompleteGammaLower(freedom / 2, value / 2);
        }
        public static double ChiSquareDistributionUpper(double freedom, double value)
        {
            return Mt.IncompleteGammaUpper(freedom / 2, value / 2);
        }

        public static double Erf(double value) { return IncompleteGammaLower(0.5, value * value) * Math.Sign(value); }


        // 不完全ベータ関数
        static double AbsMaxFPMIN(double value)
        {
            return Math.Abs(value) < DoubleFpMin ? DoubleFpMin : value;
        }
        static double betacf(double a, double b, double x)
        {
            double qab = a + b;
            double qap = a + 1;
            double qam = a - 1;
            double c = 1;
            double d = 1 / AbsMaxFPMIN(1 - qab * x / qap);
            double h = d;
            for (int m = 1; m < 10000; m++)
            {
                int m2 = 2 * m;
                double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
                d = 1 / AbsMaxFPMIN(1 + aa * d);
                c = AbsMaxFPMIN(1 + aa / c);
                h *= d * c;
                aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
                d = 1 / AbsMaxFPMIN(1 + aa * d);
                c = AbsMaxFPMIN(1 + aa / c);
                h *= d * c;
                if (Math.Abs(d * c - 1) <= DoubleEps) break;
            }
            return h;
        }

        static double betaiapprox(double a, double b, double x)
        {
            double xu;
            double a1 = a - 1;
            double b1 = b - 1;
            double mu = a / (a + b);
            double lnmu = Math.Log(mu);
            double lnmuc = Math.Log(1 - mu);
            double t = Math.Sqrt(a * b / (Mt.Sq(a + b) * (a + b + 1)));
            if (x > a / (a + b))
            {
                if (x >= 1) return 1;
                xu = Math.Min(1.0, Math.Max(mu + 10 * t, x + 5 * t));
            }
            else
            {
                if (x <= 0) return 0;
                xu = Math.Max(0.0, Math.Min(mu - 10 * t, x - 5 * t));
            }
            double sum = 0;
            for (int j = 0; j < 18; j++)
            {
                t = x + (xu - x) * Gauleg18y[j];
                sum += Gauleg18w[j] * Math.Exp(a1 * (Math.Log(t) - lnmu) + b1 * (Math.Log(1 - t) - lnmuc));
            }
            double ans = sum * (xu - x) * Math.Exp(a1 * lnmu - Mt.LogGamma(a) + b1 * lnmuc - Mt.LogGamma(b) + Mt.LogGamma(a + b));
            return ans > 0.0 ? 1.0 - ans : -ans;
        }
        public static double IncompleteBeta(double value, double param1, double param2)
        {
            if (value < 0 || value > 1) ThrowException.ArgumentOutOfRangeException("value");
            if (param1 <= 0) ThrowException.ArgumentOutOfRangeException("param1");
            if (param2 <= 0) ThrowException.ArgumentOutOfRangeException("param2");
            if (value == 0 || value == 1) return value;
            if (param1 > 3000 && param2 > 3000) return betaiapprox(param1, param2, value);
            double bt = Math.Exp(
                Mt.LogGamma(param1 + param2) - Mt.LogGamma(param1) - Mt.LogGamma(param2)
                + param1 * Math.Log(value) + param2 * Math.Log(1 - value));
            return (value < (param1 + 1) / (param1 + param2 + 2)) ?
                bt * betacf(param1, param2, value) / param1 :
                1 - bt * betacf(param2, param1, 1 - value) / param2;
        }

        // StudentT分布
        static public double StudentTPdf(double value, double freedom)
        {
            double n2 = (freedom + 1) / 2.0;
            return Mt.Gamma(n2) / (Mt.Gamma(freedom / 2.0) * Math.Sqrt(Math.PI * freedom)) * Math.Pow(1 + value * value / freedom, -n2);
        }
        // ∫∞~-t and t~∞
        static public double StudentTDistributionBilateral(double value, double freedom)
        {
            return Mt.IncompleteBeta(freedom / (freedom + value * value), freedom / 2.0, 0.5);
        }
        // ∫t~∞
        static public double StudentTDistributionUpper(double value, double freedom)
        {
            double a = Mt.StudentTDistributionBilateral(value, freedom) / 2.0;
            return value < 0.0 ? 1.0 - a : a;
        }
        // ∫-∞~t
        static public double StudentTDistributionLower(double value, double freedom)
        {
            double a = Mt.StudentTDistributionBilateral(value, freedom) / 2.0;
            return value < 0 ? a : 1 - a;
        }

        // F分布
        // ∫f~∞, Q(F|v1,v2) <0.01 で有意
        static public double FDistributionUpper(double value, double freedom1, double freedom2)
        {
            return Mt.IncompleteBeta(freedom2 / (freedom2 + freedom1 * value), freedom2 / 2.0, freedom1 / 2.0);
        }
        // ∫-∞~f = ∫0~f
        static public double FDistributionLower(double value, double freedom1, double freedom2)
        {
            return 1.0 - FDistributionUpper(value, freedom1, freedom2);
        }

        //楕円関数
        public static double EllipticTheta3(double phase, double radius)
        {
            if (radius < 0.0 || radius >= 1.0) ThrowException.ArgumentOutOfRangeException("radius");
            if (radius == 0.0) return 1.0;
            int digit = (int)Math.Ceiling(Math.Sqrt(Math.Log(DoubleEps) / Math.Log(radius)));
            double a = 0.0;
            for (int d = 1; d <= digit; d++)
                a += Math.Cos(2.0 * d * phase) * Math.Pow(radius, d * d);
            return 1.0 + 2.0 * a;
        }

        #endregion

        //---------------------------------------------------------------------------
        #region Linear functions

        // norm
        public static double Abs(IEnumerable<double> vector)
        {
            return Math.Sqrt(vector.Sum(x => x * x));
        }

        public static double Norm(IEnumerable<double> vector, double nu)
        {
            return Math.Pow(vector.Sum(x => Math.Pow(Math.Abs(x), nu)), 1 / nu);
        }
        public static double Norm1(IEnumerable<double> vector)
        {
            return vector.Sum(x => Math.Abs(x));
        }
        public static double Norm2(IEnumerable<double> vector)
        {
            return Math.Sqrt(vector.Sum(x => x * x));
        }
        public static double Norm1(Double2 vector) { return Math.Abs(vector.X) + Math.Abs(vector.Y); }
        public static double Norm1(Double3 vector) { return Math.Abs(vector.X) + Math.Abs(vector.Y) + Math.Abs(vector.Z); }
        public static double Norm2(Double2 vector) { return Math.Sqrt(vector.X * vector.X + vector.Y * vector.Y); }
        public static double Norm2(Double3 vector) { return Math.Sqrt(vector.X * vector.X + vector.Y * vector.Y + vector.Z * vector.Z); }
        public static double InnerProduct(Double2 vector) { return vector.X * vector.X + vector.Y * vector.Y; }
        public static double InnerProduct(Double3 vector) { return vector.X * vector.X + vector.Y * vector.Y + vector.Z * vector.Z; }
        public static double InnerProduct(Double2 left, Double2 right) { return left.X * right.X + left.Y * right.Y; }
        public static double InnerProduct(Double3 left, Double3 right) { return left.X * right.X + left.Y * right.Y + left.Z * right.Z; }
        public static Double2 Normalize(Double2 vector) { return vector / Norm2(vector); }
        public static Double3 Normalize(Double3 vector) { return vector / Norm2(vector); }

        public static Double3 OuterProduct(Double3 left, Double3 right)
        {
            return new Double3(
                left.Y * right.Z - left.Z * right.Y,
                left.Z * right.X - left.X * right.Z,
                left.X * right.Y - left.Y * right.X
            );
        }

        // 二次形式
        public static double InnerProduct(IList<double> vector)
        {
            double a = 0;
            for (int i = vector.Count; --i >= 0; )
                a += vector[i] * vector[i];
            return a;
        }
        public static double InnerProduct(IList<double> left, IList<double> right)
        {
            if (left.Count != right.Count) ThrowException.SizeMismatch();
            double a = 0;
            for (int i = left.Count; --i >= 0; )
                a += left[i] * right[i];
            return a;
        }

        public static double Tr(double[,] matrix)
        {
            if (matrix.GetLength(0) != matrix.GetLength(1)) ThrowException.SizeMismatch();
            double a = 0;
            for (int i = matrix.GetLength(0); --i >= 0; )
                a += matrix[i, i];
            return a;
        }

        public static double QuadraticForm(double[,] matrix, IList<double> vector)
        {
            if (matrix.GetLength(0) != vector.Count || matrix.GetLength(1) != vector.Count) ThrowException.SizeMismatch();
            double a = 0;
            for (int i = matrix.GetLength(0); --i >= 0; )
            {
                double b = 0;
                for (int j = matrix.GetLength(1); --j >= 0; )
                    b += matrix[i, j] * vector[j];
                a += b * vector[i];
            }
            return a;
        }
        public static double QuadraticForm(IList<double> left, double[,] matrix, IList<double> right)
        {
            if (matrix.GetLength(0) != left.Count || matrix.GetLength(1) != right.Count) ThrowException.SizeMismatch();
            double a = 0;
            for (int i = matrix.GetLength(0); --i >= 0; )
            {
                double b = 0;
                for (int j = matrix.GetLength(1); --j >= 0; )
                    b += matrix[i, j] * right[j];
                a += b * left[i];
            }
            return a;
        }
        public static double TrMul(double[,] matrix1, double[,] matrix2)
        {
            if (matrix1.GetLength(1) != matrix2.GetLength(0) || matrix1.GetLength(0) != matrix2.GetLength(1)) ThrowException.SizeMismatch();
            double a = 0;
            for (int i = matrix1.GetLength(0); --i >= 0; )
                for (int j = matrix1.GetLength(1); --j >= 0; )
                    a += matrix1[i, j] * matrix2[j, i];
            return a;
        }
        public static double TrMul(double[,] matrix1, double[,] matrix2, double[,] matrix3)
        {
            if (matrix1.GetLength(1) != matrix2.GetLength(0) || matrix2.GetLength(1) != matrix3.GetLength(0) || matrix1.GetLength(0) != matrix3.GetLength(1)) ThrowException.SizeMismatch();
            double a = 0;
            for (int i = matrix2.GetLength(0); --i >= 0; )
                for (int j = matrix2.GetLength(1); --j >= 0; )
                {
                    double b = 0;
                    for (int k = matrix1.GetLength(0); --k >= 0; )
                        b += matrix1[k, i] * matrix3[j, k];
                    a += b * matrix2[i, j];
                }
            return a;
        }

        // 行列
        public static double[,] UnitMatrix(int size)
        {
            double[,] result = new double[size, size];
            for (int i = size; --i >= 0; )
                result[i, i] = 1;
            return result;
        }
        public static double[,] DiagonalMatrix(int size, double value)
        {
            double[,] result = new double[size, size];
            for (int i = size; --i >= 0; )
                result[i, i] = value;
            return result;
        }
        public static double[,] DiagonalMatrix(IList<double> values)
        {
            double[,] result = new double[values.Count, values.Count];
            for (int i = values.Count; --i >= 0; )
                result[i, i] = values[i];
            return result;
        }
        public static double[,] VectorsToMatrix(IList<double> vector)
        {
            double[,] result = new double[vector.Count, vector.Count];
            for (int i = vector.Count; --i >= 0; )
                for (int j = vector.Count; --j >= 0; )
                    result[i, j] = vector[i] * vector[j];
            return result;
        }
        public static double[,] VectorsToMatrix(IList<double> left, IList<double> right)
        {
            double[,] result = new double[left.Count, right.Count];
            for (int i = 0; i < left.Count; i++)
                for (int j = 0; j < right.Count; j++)
                    result[i, j] = left[i] * right[j];
            return result;
        }
        public static void LetSymmetrical(double[,] matrix)
        {
            if (matrix.GetLength(0) != matrix.GetLength(1)) ThrowException.SizeMismatch();
            for (int i = matrix.GetLength(0); --i >= 0; )
                for (int j = i; --j >= 0; )
                {
                    double a = (matrix[i, j] + matrix[j, i]) * 0.5;
                    matrix[i, j] = a;
                    matrix[j, i] = a;
                }
        }

        // 行列計算
        public static double[] Neg(IList<double> vector)
        {
            double[] result = new double[vector.Count];
            for (int i = result.Length; --i >= 0; )
                result[i] = -vector[i];
            return result;
        }
        public static void LetNeg(IList<double> vector)
        {
            for (int i = vector.Count; --i >= 0; )
                vector[i] *= -1;
        }
        public static double[] Add(IList<double> left, IList<double> right)
        {
            if (left.Count != right.Count) ThrowException.SizeMismatch();
            double[] result = new double[left.Count];
            for (int i = result.Length; --i >= 0; )
                result[i] = left[i] + right[i];
            return result;
        }
        public static void LetAdd(IList<double> left, IList<double> right)
        {
            if (left.Count != right.Count) ThrowException.SizeMismatch();
            for (int i = left.Count; --i >= 0; )
                left[i] += right[i];
        }
        public static double[] Sub(IList<double> left, IList<double> right)
        {
            if (left.Count != right.Count) ThrowException.SizeMismatch();
            double[] result = new double[left.Count];
            for (int i = result.Length; --i >= 0; )
                result[i] = left[i] - right[i];
            return result;
        }
        public static void LetSub(IList<double> left, IList<double> right)
        {
            if (left.Count != right.Count) ThrowException.SizeMismatch();
            for (int i = left.Count; --i >= 0; )
                left[i] -= right[i];
        }

        public static double[,] Neg(this double[,] matrix)
        {
            double[,] result = new double[matrix.GetLength(0), matrix.GetLength(1)];
            for (int i = result.GetLength(0); --i >= 0; )
                for (int j = result.GetLength(1); --j >= 0; )
                    result[i, j] = -matrix[i, j];
            return result;
        }
        public static void LetNeg(this double[,] matrix)
        {
            for (int i = matrix.GetLength(0); --i >= 0; )
                for (int j = matrix.GetLength(1); --j >= 0; )
                    matrix[i, j] *= -1;
        }
        public static double[,] Add(double[,] left, double[,] right)
        {
            if (left.GetLength(0) != right.GetLength(0) || left.GetLength(1) != right.GetLength(1)) ThrowException.SizeMismatch();
            double[,] result = new double[left.GetLength(0), left.GetLength(1)];
            for (int i = result.GetLength(0); --i >= 0; )
                for (int j = result.GetLength(1); --j >= 0; )
                    result[i, j] = left[i, j] + right[i, j];
            return result;
        }
        public static void LetAdd(double[,] left, double[,] right)
        {
            if (left.GetLength(0) != right.GetLength(0) || left.GetLength(1) != right.GetLength(1)) ThrowException.SizeMismatch();
            for (int i = left.GetLength(0); --i >= 0; )
                for (int j = left.GetLength(1); --j >= 0; )
                    left[i, j] += right[i, j];
        }
        public static double[,] Sub(double[,] left, double[,] right)
        {
            if (left.GetLength(0) != right.GetLength(0) || left.GetLength(1) != right.GetLength(1)) ThrowException.SizeMismatch();
            double[,] result = new double[left.GetLength(0), left.GetLength(1)];
            for (int i = result.GetLength(0); --i >= 0; )
                for (int j = result.GetLength(1); --j >= 0; )
                    result[i, j] = left[i, j] - right[i, j];
            return result;
        }
        public static void LetSub(double[,] left, double[,] right)
        {
            if (left.GetLength(0) != right.GetLength(0) || left.GetLength(1) != right.GetLength(1)) ThrowException.SizeMismatch();
            for (int i = left.GetLength(0); --i >= 0; )
                for (int j = left.GetLength(1); --j >= 0; )
                    left[i, j] -= right[i, j];
        }

        public static double[] Mul(double left, IList<double> right) { return Mul(right, left); }
        public static double[] Mul(IList<double> left, double right)
        {
            double[] result = new double[left.Count];
            for (int i = result.Length; --i >= 0; )
                result[i] = left[i] * right;
            return result;
        }
        public static void LetMul(IList<double> left, double right)
        {
            for (int i = left.Count; --i >= 0; )
                left[i] *= right;
        }

        public static double[,] Mul(double left, double[,] right) { return Mul(right, left); }
        public static double[,] Mul(double[,] left, double right)
        {
            double[,] result = new double[left.GetLength(0), left.GetLength(1)];
            for (int i = result.GetLength(0); --i >= 0; )
                for (int j = result.GetLength(1); --j >= 0; )
                    result[i, j] = left[i, j] * right;
            return result;
        }
        public static void LetMul(double[,] left, double right)
        {
            for (int i = left.GetLength(0); --i >= 0; )
                for (int j = left.GetLength(1); --j >= 0; )
                    left[i, j] *= right;
        }

        public static double[] Mul(double[,] left, IList<double> right)
        {
            if (left.GetLength(1) != right.Count) ThrowException.SizeMismatch();
            double[] result = new double[left.GetLength(0)];
            for (int i = result.Length; --i >= 0; )
            {
                double a = 0;
                for (int j = left.GetLength(1); --j >= 0; )
                    a += left[i, j] * right[j];
                result[i] = a;
            }
            return result;
        }
        public static double[] Mul(IList<double> left, double[,] right)
        {
            if (left.Count != right.GetLength(0)) ThrowException.SizeMismatch();
            double[] result = new double[right.GetLength(1)];
            for (int j = left.Count; --j >= 0; )
            {
                double val = left[j];
                for (int i = result.Length; --i >= 0; )
                    result[i] += val * right[j, i];
            }
            return result;
        }

        public static double[,] Mul(double[,] left, double[,] right)
        {
            if (left.GetLength(1) != right.GetLength(0)) ThrowException.SizeMismatch();
            double[,] result = new double[left.GetLength(0), right.GetLength(1)];
            for (int i = result.GetLength(0); --i >= 0; )
                for (int k = left.GetLength(1); --k >= 0; )
                {
                    double val = left[i, k];
                    for (int j = result.GetLength(1); --j >= 0; )
                        result[i, j] += val * right[k, j];
                }
            return result;
        }
        //↓のunsafe関数を使いたくない場合は↑の関数を使う
        //public unsafe static double[,] Mul(double[,] left, double[,] right)
        //{
        //    int I = left.GetLength(0);
        //    int J = right.GetLength(1);
        //    int K = right.GetLength(0);
        //    if (K != left.GetLength(1)) ThrowException.SizeMismatch();
        //    double[,] result = new double[I, J];
        //    int U = J & ~3;
        //    fixed (double* qr = &right[0, 0])
        //    fixed (double* qz = &result[0, 0])
        //    {
        //        for (int i = I; --i >= 0; )
        //        {
        //            double* pz = qz + i * J;
        //            for (int k = K; --k >= 0; )
        //            {
        //                double* pr = qr + k * J;
        //                double val = left[i, k];
        //                for (int j = J; --j >= U; )
        //                    *(pz + j) += val * *(pr + j);
        //                for (int j = U; (j -= 4) >= 0; )
        //                {
        //                    *(pz + j + 3) += val * *(pr + j + 3);
        //                    *(pz + j + 2) += val * *(pr + j + 2);
        //                    *(pz + j + 1) += val * *(pr + j + 1);
        //                    *(pz + j + 0) += val * *(pr + j + 0);
        //                }
        //            }
        //        }
        //    }
        //    return result;
        //}

        public static double[,] T(this double[,] matrix)
        {
            double[,] result = new double[matrix.GetLength(1), matrix.GetLength(0)];
            for (int i = result.GetLength(0); --i >= 0; )
                for (int j = result.GetLength(1); --j >= 0; )
                    result[i, j] = matrix[j, i];
            return result;
        }

        static double Det1by1(double[,] matrix)
        {
            return matrix[0, 0];
        }
        static double[,] Inverse1by1(double[,] matrix)
        {
            double det = matrix[0, 0];
            if (det == 0) ThrowException.WriteLine("Inverse: singular matrix");
            return new double[1, 1]{
                {1.0 / det}
            };
        }
        static double Det2by2(double[,] matrix)
        {
            double a = matrix[0, 0], b = matrix[0, 1];
            double c = matrix[1, 0], d = matrix[1, 1];
            return a * d - b * c;
        }
        static double[,] Inverse2by2(double[,] matrix)
        {
            double a = matrix[0, 0], b = matrix[0, 1];
            double c = matrix[1, 0], d = matrix[1, 1];
            double det = a * d - b * c;
            if (det == 0) ThrowException.WriteLine("Inverse: singular matrix");
            return new double[2, 2]{
                {d / det, -b / det},
                {-c / det, a / det}
            };
        }
        static double Det3by3(double[,] matrix)
        {
            double a = matrix[0, 0], b = matrix[0, 1], c = matrix[0, 2];
            double d = matrix[1, 0], e = matrix[1, 1], f = matrix[1, 2];
            double g = matrix[2, 0], h = matrix[2, 1], i = matrix[2, 2];
            return a * (e * i - f * h) + b * (f * g - d * i) + c * (d * h - e * g);
        }
        static double[,] Inverse3by3(double[,] matrix)
        {
            double a = matrix[0, 0], b = matrix[0, 1], c = matrix[0, 2];
            double d = matrix[1, 0], e = matrix[1, 1], f = matrix[1, 2];
            double g = matrix[2, 0], h = matrix[2, 1], i = matrix[2, 2];
            double eifh = e * i - f * h, fgdi = f * g - d * i, dheg = d * h - e * g;
            double det = a * eifh + b * fgdi + c * dheg;
            if (det == 0) ThrowException.WriteLine("Inverse: singular matrix");
            return new double[3, 3]{
                { eifh / det, (c * h - b * i) / det, (b * f - c * e) / det},
                { fgdi / det, (a * i - c * g) / det, (c * d - a * f) / det},
                { dheg / det, (b * g - a * h) / det, (a * e - b * d) / det}
            };
        }

        static int[] LUDecomposition(double[][] matrix)
        {
            int n = matrix.GetLength(0);
            int[] pivot = new int[n];
            double[] vv = new double[n];
            for (int i = 0; i < n; i++)
            {
                double max = Ex.Range(n).Select(j => Math.Abs(matrix[i][j])).Max();
                if (max == 0.0) ThrowException.WriteLine("LUDecomposition: singular matrix");
                vv[i] = 1.0 / max;
            }

            for (int k = 0; k < n; k++)
            {
                int p = k + Ex.FromTo(k, n).Select(i => vv[i] * Math.Abs(matrix[i][k])).MaxIndex();
                if (p != k)
                {
                    Mt.Swap(ref matrix[p], ref matrix[k]);
                    vv[p] = vv[k];
                }
                pivot[k] = p;
                var mk = matrix[k];
                if (mk[k] == 0.0) { ThrowException.WriteLine("LUDecomposition: singular matrix"); mk[k] = 1.0e-40; }
                for (int i = k; ++i < n; )  //順不同並列可
                {
                    var mi = matrix[i];
                    double temp = mi[k] /= mk[k];
                    for (int j = k; ++j < n; )
                        mi[j] -= temp * mk[j];
                }
            }
            return pivot;
        }

        public static double[,] Inverse(this double[,] matrix)
        {
            int n = matrix.GetLength(0);
            if (n == 1) return Inverse1by1(matrix);
            if (n == 2) return Inverse2by2(matrix);
            if (n == 3) return Inverse3by3(matrix);

            double[][] LU = New.Array(n, i => New.Array(n, j => matrix[j, i]));
            int[] pivot = LUDecomposition(LU);
            int[] index = New.Array(n, i => i);
            for (int i = 0; i < n; i++)
                Mt.Swap(ref index[i], ref index[pivot[i]]);

            double[,] result = new double[n, n];
            double[] vec = new double[n];
            for (int i = n; --i >= 0; )  //順不同並列可
            {
                vec.Clear();
                vec[i] = 1.0;
                for (int j = i; ++j < n; )
                {
                    double sum = 0;
                    for (int k = i; k < j; k++)
                        sum -= LU[j][k] * vec[k];
                    vec[j] = sum;
                }
                for (int j = n; --j >= 0; )
                {
                    double sum = vec[j];
                    for (int k = j; ++k < n; )
                        sum -= LU[j][k] * vec[k];
                    vec[j] = sum / LU[j][j];
                }
                int idx = index[i];
                for (int j = n; --j >= 0; )
                    result[idx, j] = vec[j];
            }
            return result;
        }

        public static double Det(this double[,] matrix)
        {
            int n = matrix.GetLength(0);
            if (n == 1) return Det1by1(matrix);
            if (n == 2) return Det2by2(matrix);
            if (n == 3) return Det3by3(matrix);

            double[][] L = New.Array(n, i => New.Array(n, j => matrix[i, j]));
            int[] pivot = LUDecomposition(L);
            int swap = Ex.Range(n).Count(i => pivot[i] != i);
            double sign = swap % 2 == 0 ? 1.0 : -1.0;
            return sign * Ex.Range(n).Product(i => L[i][i]);
        }

        // Inverse of symmetric positive definite matrix
        static double[,] CholeskyDecomposition(double[,] matrix)
        {
            int n = matrix.GetLength(0);
            double[,] L = new double[n, n];
            //対角と下三角部分を計算
            for (int j = 0; j < n; j++)
            {
                for (int i = j; i < n; i++)
                {
                    double sum = matrix[j, i];
                    for (int k = j; --k >= 0; )
                        sum -= L[j, k] * L[i, k];
                    if (j == i && sum <= 0.0) ThrowException.WriteLine("CholeskyDecomposition: non-positive definite matrix");
                    L[i, j] = (j == i) ? Math.Sqrt(sum) : sum / L[j, j];
                }
            }
            return L;
        }
        public static double[,] InverseSymmetricPositiveDefinite(double[,] matrix)
        {
            int n = matrix.GetLength(0);
            double[,] L = CholeskyDecomposition(matrix);
            //対角と上三角部分を計算
            for (int j = 0; j < n; j++)
            {
                L[j, j] = 1.0 / L[j, j];
                for (int i = j; --i >= 0; )
                {
                    double sum = 0.0;
                    for (int k = j; --k >= i; )
                        sum -= L[j, k] * L[i, k];
                    L[i, j] = sum * L[j, j];
                }
            }
            for (int j = n; --j >= 0; )
            {
                for (int i = 0; i <= j; i++)
                {
                    double sum = L[i, j];
                    for (int k = j; ++k < n; )
                        sum -= L[k, j] * L[i, k];
                    L[i, j] = sum * L[j, j];
                }
            }
            //下三角部分を上三角部分からコピー
            for (int i = n; --i >= 0; )
                for (int j = i; --j >= 0; )
                    L[i, j] = L[j, i];
            return L;
        }


        //Sqrt(x*x + y*y) をoverflow, underflowに気をつけて計算
        static double pythag(double x, double y)
        {
            double z = Math.Sqrt(x * x + y * y);
            if (!double.IsPositiveInfinity(z)) return z;
            return Math.Abs(x) < Math.Abs(y) ?
                Math.Sqrt(1 + Mt.Sq(x / y)) * Math.Abs(y) :
                Math.Sqrt(1 + Mt.Sq(y / x)) * Math.Abs(x);
        }
        static void _QR(double[] vector0, double[] vector1, double c, double s)
        {
            for (int i = vector0.Length; --i >= 0; )
            {
                double x = vector0[i];
                double y = vector1[i];
                vector0[i] = x * c + y * s;
                vector1[i] = y * c - x * s;
            }
        }

        //M: symmetric matrix
        //Mを転置させた計算の方が速いため改変
        static double[] Householder(double[][] M, double[] D)
        {
            int N = M.Length;
            double[] E = new double[N];
            for (int i = N; --i >= 2; )
            {
                double[] Mi = new double[i];
                for (int l = i; --l >= 0; ) Mi[l] = M[l][i];
                double scale = Mi.Sum(x => Math.Abs(x));
                if (scale == 0) continue;
                for (int l = i; --l >= 0; ) Mi[l] /= scale;

                double Di;
                {
                    double f = Mi.Sum(x => x * x);
                    double g = Math.Sqrt(f);
                    if (Mi[i - 1] < 0) g *= -1;
                    E[i] = -g * scale;
                    D[i] = Di = f + g * Mi[i - 1];
                    Mi[i - 1] += g;
                }

                double hh = 0;
                for (int j = 0; j < i; j++)
                {
                    M[i][j] = Mi[j] / Di;
                    {
                        double a = 0;
                        for (int l = i; --l > j; ) a += M[j][l] * Mi[l];
                        for (int l = j + 1; --l >= 0; ) a += M[l][j] * Mi[l];
                        E[j] = a / Di;
                    }
                    hh += E[j] * Mi[j];
                }
                hh /= 2 * Di;
                for (int j = 0; j < i; j++)
                {
                    E[j] -= Mi[j] * hh;
                    for (int k = j + 1; --k >= 0; ) M[k][j] -= Mi[j] * E[k] + Mi[k] * E[j];
                }
                for (int l = i; --l >= 0; ) M[l][i] = Mi[l];
            }
            E[1] = M[0][1];

            for (int i = 0; i < N; i++)
            {
                if (D[i] != 0)
                {
                    for (int j = 0; j < i; j++)
                    {
                        double a = 0;
                        for (int l = i; --l >= 0; ) a += M[j][l] * M[l][i];
                        for (int l = i; --l >= 0; ) M[j][l] -= M[i][l] * a;
                    }
                }
                D[i] = M[i][i];
                M[i][i] = 1;
                for (int l = i; --l >= 0; ) M[i][l] = 0;
                for (int l = i; --l >= 0; ) M[l][i] = 0;
            }
            return E;
        }
        static void QLMethod(double[][] M, double[] D, double[] E)
        {
            int N = M.Length;
            for (int l = 0; l < N; l++)
            {
                while (true)
                {
                    int m;
                    for (m = l; m < N - 1; m++)
                        if (IsTooSmall(E[m], Math.Abs(D[m]) + Math.Abs(D[m + 1])) || Math.Abs(E[m]) < Mt.DoubleEpsilon) break;
                    if (m == l) break;

                    double g = (D[l + 1] - D[l]) / E[l] * 0.5;
                    g += g * Math.Sqrt(1 + 1 / (g * g));
                    g = E[l] / g + D[m] - D[l];
                    double s = 1, c = 1, p = 0;
                    int i;
                    for (i = m; --i >= l; )
                    {
                        double b = E[i] * c;
                        {
                            double f = E[i] * s;
                            double py = pythag(f, g);
                            E[i + 1] = py;
                            if (py == 0) break;
                            s = f / py;
                            c = g / py;
                        }
                        {
                            double q = D[i + 1] - p;
                            double r = (D[i] - q) * s + b * c * 2;
                            p = r * s;
                            D[i + 1] = q + p;
                            g = r * c - b;
                        }
                        _QR(M[i], M[i + 1], c, -s);
                    }
                    if (i < l) E[l] = g;
                    D[i + 1] -= p;
                    E[m] = 0;
                }
            }
        }
        //matrix: symmetric matrix
        //matrix => V * diag[D] * V^T
        //固有値の配列と固有ベクトルの配列を返す
        public static Tuple<double[], double[][]> EigenValueDecomposition(double[,] matrix)
        {
            int N = matrix.GetLength(0);
            var D = new double[N];
            var V = New.Array(N, i => New.Array(N, j => matrix[i, j]));

            double[] E = Householder(V, D);
            for (int i = 0; i < N - 1; i++) E[i] = E[i + 1];
            E[N - 1] = 0;
            QLMethod(V, D, E);

            int[] order = D.SortIndex(false);
            D.SortAsIndex(order);
            V.SortAsIndex(order);
            for (int i = N; --i >= 0; )
                if (V[i].Count(a => a < 0) > N / 2) Mt.LetMul(V[i], -1);
            return new Tuple<double[], double[][]>(D, V);
        }
        //matrix => V * diag[D] * V^T
        //V, Dを返す
        public static Tuple<double[,], double[]> EigenValueDecompositionM(double[,] matrix)
        {
            var T = EigenValueDecomposition(matrix);
            return new Tuple<double[,], double[]>(
                New.Array(T.Item2[0].Length, T.Item2.Length, (i, j) => T.Item2[j][i]),
                T.Item1
            );
        }

        //matrix -> U diag[W] V^T
        //U, V: 縦ベクトルの配列を返す
        public static Tuple<double[][], double[], double[][]> SingularValueDecomposition(double[,] matrix)
        {
            int XN = matrix.GetLength(1);
            int YN = matrix.GetLength(0);
            int ZN = Math.Min(YN, XN);
            var U = New.Array(XN, i => New.Array(YN, j => matrix[j, i]));
            var V = New.Array(XN, i => new double[XN]);
            var W = new double[XN];
            var R = new double[XN];
            double anorm = 0;
            for (int i = 0; i < XN; i++)
            {
                var I = Ex.FromTo(i, YN);
                if (i < YN)
                {
                    double scale = I.Sum(k => Math.Abs(U[i][k]));
                    if (scale != 0.0)
                    {
                        I.ForEach(k => { U[i][k] /= scale; });
                        double s = I.Sum(k => Mt.Sq(U[i][k]));
                        double f = U[i][i];
                        double g = Math.Sqrt(s);
                        if (f > 0) g *= -1;
                        double h = f * g - s;
                        U[i][i] = f - g;
                        for (int j = i; ++j < XN; )
                        {
                            double coef = 0;
                            for (int k = i; k < YN; k++) coef += U[j][k] * U[i][k];
                            coef /= h;
                            for (int k = i; k < YN; k++) U[j][k] += coef * U[i][k];
                        }
                        I.ForEach(k => { U[i][k] *= scale; });
                        W[i] = scale * g;
                    }
                }
                int l = i + 1;
                var L = Ex.FromTo(l, XN);
                if (i < YN && i < XN - 1)
                {
                    double scale = L.Sum(k => Math.Abs(U[k][i]));
                    if (scale != 0.0)
                    {
                        L.ForEach(k => { U[k][i] /= scale; });
                        double s = L.Sum(k => Mt.Sq(U[k][i]));
                        double f = U[l][i];
                        double g = -Math.Sign(f) * Math.Sqrt(s);
                        double h = f * g - s;
                        U[l][i] = f - g;
                        L.ForEach(k => { R[k] = U[k][i] / h; });
                        for (int j = l; j < YN; j++)
                        {
                            double coef = 0;
                            for (int k = l; k < XN; k++) coef += U[k][j] * U[k][i];
                            for (int k = l; k < XN; k++) U[k][j] += coef * R[k];
                        }
                        L.ForEach(k => { U[k][i] *= scale; });
                        R[l] = scale * g;
                    }
                }
                Mt.LetMax(ref anorm, Math.Abs(W[i]) + Math.Abs(R[i]));
            }

            for (int i = XN; --i >= 0; )
            {
                if (i < XN - 1)
                {
                    int l = i + 1;
                    var L = Ex.FromTo(l, XN);
                    double g = R[l];
                    if (g != 0.0)
                    {
                        if (U[l][i] != 0.0)
                        {
                            double coef = 1 / U[l][i] / g;
                            L.ForEach(k => { V[i][k] = U[k][i] * coef; });
                        }
                        for (int j = l; j < XN; j++)
                        {
                            double coef = 0;
                            for (int k = l; k < XN; k++) coef += V[j][k] * U[k][i];
                            for (int k = l; k < XN; k++) V[j][k] += coef * V[i][k];
                        }
                    }
                    L.ForEach(k => { V[i][k] = 0.0; });
                    L.ForEach(k => { V[k][i] = 0.0; });
                }
                V[i][i] = 1.0;
            }

            for (int i = ZN; --i >= 0; )
            {
                int l = i + 1;
                var I = Ex.FromTo(i, YN);
                Ex.FromTo(l, XN).ForEach(k => { U[k][i] = 0.0; });
                double g = W[i];
                if (g != 0.0)
                {
                    g = 1 / g;
                    for (int j = l; j < XN; j++)
                    {
                        double coef = 0;
                        for (int k = l; k < YN; k++) coef += U[j][k] * U[i][k];
                        coef = coef / U[i][i] * g;
                        for (int k = i; k < YN; k++) U[j][k] += coef * U[i][k];  //
                    }
                    I.ForEach(k => { U[i][k] *= g; });
                }
                else
                    I.ForEach(k => { U[i][k] = 0.0; });
                U[i][i] += 1;
            }

            for (int k = XN; --k >= 0; )
            {
                for (int its = 30; ; its--)
                {
                    int l;
                    for (l = k; l > 0; l--)
                    {  // R[0]==0;
                        if (IsTooSmall(R[l], anorm)) break;
                        if (IsTooSmall(W[l - 1], anorm))
                        {
                            double c = 0, s = 1;
                            for (int i = l; i <= k; i++)
                            {
                                double f = s * R[i];
                                R[i] *= c;
                                if (IsTooSmall(f, anorm)) break;
                                double g = W[i];
                                double py = pythag(f, g);
                                W[i] = py;
                                c = g / py;
                                s = -f / py;
                                _QR(U[l - 1], U[i], c, s);
                            }
                            break;
                        }
                    }

                    if (l == k)
                    {
                        if (W[k] < 0) { W[k] *= -1; Mt.LetMul(V[k], -1); }
                        break;
                    }
                    if (its == 0) { Console.WriteLine("DecompositionSingularValue: iteration>30"); break; }

                    {
                        double f = 0.0;
                        double x = W[l];
                        {
                            double y = W[k - 1];
                            double z = W[k];
                            double g = R[k - 1];
                            double h = R[k];
                            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2 * h * y);
                            double py = f * (1 + Math.Sqrt(1 + 1 / (f * f)));
                            f = ((x - z) * (x + z) + h * (y / py - h)) / x;
                        }
                        double c = 1, s = 1;
                        for (int j = l; j < k; j++)
                        {
                            int i = j + 1;
                            double g = R[i] * c;
                            double h = R[i] * s;
                            {
                                double py = pythag(f, h);
                                R[j] = py;
                                c = f / py;
                                s = h / py;
                            }
                            f = x * c + g * s;
                            g = g * c - x * s;

                            double y = W[i] * c;
                            h = W[i] * s;
                            _QR(V[j], V[i], c, s);
                            {
                                double py = pythag(f, h);
                                W[j] = py;
                                if (py != 0.0)
                                {
                                    c = f / py;
                                    s = h / py;
                                }
                            }
                            f = g * c + y * s;
                            x = y * c - g * s;
                            _QR(U[j], U[i], c, s);
                        }
                        R[l] = 0;
                        R[k] = f;
                        W[k] = x;
                    }
                }
            }
            {
                var order = W.SortIndex((x, y) => x < y ? +1 : -1);
                W.SortAsIndex(order);
                U.SortAsIndex(order);
                V.SortAsIndex(order);

                for (int k = 0; k < XN; k++)
                {
                    int s = U[k].Count(a => a < 0) + V[k].Count(a => a < 0);
                    if (s > (YN + XN) / 2)
                    {
                        Mt.LetMul(U[k], -1);
                        Mt.LetMul(V[k], -1);
                    }
                }
            }
            return new Tuple<double[][], double[], double[][]>(U, W, V);
        }
        //matrix -> U diag[W] V^T
        //U, W, Vを返す
        public static Tuple<double[,], double[], double[,]> SingularValueDecompositionM(double[,] matrix)
        {
            var T = SingularValueDecomposition(matrix);
            return new Tuple<double[,], double[], double[,]>(
                New.Array(T.Item1[0].Length, T.Item1.Length, (i, j) => T.Item1[j][i]),
                T.Item2,
                New.Array(T.Item3.Length, T.Item3.Length, (i, j) => T.Item3[j][i])
            );
        }


        // Rotation Matrix
        public static Double2 Rotate(Double2 vector, double angle)
        {
            double cos = Math.Cos(angle);
            double sin = Math.Sin(angle);
            return new Double2(cos * vector.X - sin * vector.Y, sin * vector.X + cos * vector.Y);
        }
        // 右手座標系ならvectorをaxisを軸にして右ねじの回転方向に回転させる．axisの長さ=1
        public static Double3 Rotate(Double3 vector, Double3 axis, double angle)
        {
            double cos = Math.Cos(angle);
            double sin = Math.Sin(angle);
            double cos1 = 1 - cos;
            double x = axis.X;
            double y = axis.Y;
            double z = axis.Z;
            return new Double3(
                (cos1 * x * x + cos) * vector.X + (cos1 * x * y - sin * z) * vector.Y + (cos1 * z * x + sin * y) * vector.Z,
                (cos1 * x * y + sin * z) * vector.X + (cos1 * y * y + cos) * vector.Y + (cos1 * y * z - sin * x) * vector.Z,
                (cos1 * z * x - sin * y) * vector.X + (cos1 * y * z + sin * x) * vector.Y + (cos1 * z * z + cos) * vector.Z
            );
        }
        public static double[,] RotationMatrix(double angle)
        {
            double cos = Math.Cos(angle);
            double sin = Math.Sin(angle);
            double[,] matrix = new double[2, 2]{
                { cos, -sin },
                { sin, cos }
            };
            return matrix;
        }
        public static double[,] RotationMatrix(Double3 axis, double angle)
        {
            double cos = Math.Cos(angle);
            double sin = Math.Sin(angle);
            double cos1 = 1 - cos;
            double x = axis.X;
            double y = axis.Y;
            double z = axis.Z;
            double[,] matrix = new double[3, 3] {
                { cos1 * x * x + cos, cos1 * x * y - sin * z, cos1 * z * x + sin * y },
                { cos1 * x * y + sin * z, cos1 * y * y + cos, cos1 * y * z - sin * x },
                { cos1 * z * x - sin * y, cos1 * y * z + sin * x, cos1 * z * z + cos },
            };
            return matrix;
        }

        #endregion

    }

}

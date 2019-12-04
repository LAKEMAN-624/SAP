using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;

namespace dealcharacter
{
    class Program
    {
        static void Main(string[] args)
        {
            // C#中的字符串分割，并转存为ASCII码数组

            String c = "abcdefg";

            ASCIIEncoding ascii = new ASCIIEncoding();

            int[] num = new int[c.Length];

            for (int i = 0; i < c.Length; i++)

            {

                num[i] = (int)ascii.GetBytes(c)[i]; //通过ASCIIEncoding类的对象调用GetBytes方法将字符串转变成字符，然后将字符转变成8位的Byte类型，再将Byte类型转变成int类型。

            }

            int m = c.Length;

            FortranMethod.S(ref num[0], ref m);//在C#中使用ref关键字，表示参数传递时使用引用类型
        }
    }
    public static class FortranMethod
    {
        [DllImport("character.dll", SetLastError = true, CharSet = CharSet.Unicode,CallingConvention = CallingConvention.Cdecl)]//DLL程序在C#代码中的入口点

        public static extern void S(ref int ka, ref int m);
    }
}

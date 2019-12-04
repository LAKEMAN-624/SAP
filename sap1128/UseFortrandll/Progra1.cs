using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;

namespace UseFortrandll
{
    
    class Program1
    {


        static void Main(string[] args)
        {

            Console.WriteLine("请输入两个数相加：");
            double num1 = Convert.ToDouble(Console.ReadLine());

            double num2 = Convert.ToDouble(Console.ReadLine());
            Console.WriteLine("输入的两个数是：" + num1 + " ," + num2);
            double sum = FortranMethod.Add(num1, num2);
            Console.WriteLine("求和结果是：" + sum);

            double[] Array = { 1, 5, 2, 4, 3, 7, 6 };
            Console.WriteLine("初始数组：");
            for (int i = 0; i < Array.Length; i++)
                Console.Write(Array[i] + " ");

            double b = FortranMethod.Sortandfindmax(Array, Array.Length);
            Console.WriteLine("\n" + "排序后：");

            for (int i = 0; i < Array.Length; i++)
                Console.Write(Array[i] + " ");

            Console.WriteLine("\n" + "最大值为：");
            Console.WriteLine(b);

            Console.ReadKey();

        }
      
    }
    
}
public static class FortranMethod
{
    [DllImport("ArrayDll.dll", CallingConvention = CallingConvention.Cdecl)]
    public static extern double Add(double a, double b);

    [DllImport("ArrayDll.dll", CallingConvention = CallingConvention.Cdecl)]
    public static extern double Sortandfindmax(double[] array, int length);
}
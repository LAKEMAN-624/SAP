using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;
       
namespace UseFortrandll
{

    class Program
    {
       
         

        static void Main(string[] args)
        {            
            double[] angle =new double[] { 30.0 ,45.0};
            int speed = 2;
            double[] distance = new double[2];
          
            FortranMethod.ADDRESULT (ref angle, ref speed, ref distance);
           Console.WriteLine(distance);
           
            Console.ReadKey();           
        }       
    }
    public static class FortranMethod
    {
         [DllImport("throw.dll", SetLastError = true, CharSet = CharSet.Unicode, CallingConvention = CallingConvention.Cdecl)] 
        public static extern double ADDRESULT(ref double[] A, ref int B, ref double[] C);
    }
}

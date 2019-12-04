using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;
namespace DEALRPRM
{
    class Program
    {
        static void Main(string[] args)
        {
            int IA = 800;
            int IC = 200;
            double[] BB = new double[200];
            int NRWK = 20;
            double[] RPRM = new double[10] { 1,2,3,4,5,6,7,8,9,10};
            FortranMethod.DOOT(ref  IA, ref RPRM[0], ref NRWK, ref  BB[0], ref  IC);

        }
    }
    public static class FortranMethod
    {
        [DllImport("RPRM.dll", CallingConvention = CallingConvention.StdCall)]
        //[DllImport("DOOT1.dll")]
        //public static extern double DOT(ref int INFO, ref int METHOD, ref int IPRINT, ref int NDV, ref int NCON, ref double[] X, ref double[] XL, ref double[] XU, ref double OBJ, ref int MINMAX, ref double[] G, ref double[] RPRM, ref int[] IPRM, ref double[] WK, ref int NRWK, ref int[] IWK, ref int NRIWK);
        //public static extern double DOT( int INFO,  int METHOD,  int IPRINT,  int NDV,  int NCON, double[] X,  double[] XL, double[] XU,  double OBJ,  int MINMAX,  double[] G, double[] RPRM, int[] IPRM,  double[] WK, int NRWK,  int[] IWK, int NRIWK);
        public static extern double DOOT(ref int IA, ref double RPRM, ref int NRWK, ref double BB, ref int IC);

       
    }
  
}

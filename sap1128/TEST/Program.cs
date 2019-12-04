using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;
namespace TEST
{
    class Program
    {
        static void Main(string[] args)
        {
            double[] X = new double[2]; double[] XL = new double[2];
            double[] XU = new double[2]; double[] G = new double[2];
            double[] WK = new double[800];
            int[] IWK = new int[200];
            double[] RPRM = new double[20];
            int[] IPRM = new int[20];
            int NRWK = 800; int NRIWK = 200;
            double OBJ = new double();
            OBJ = 0;
            for (int i = 0; i < 20; i++)
            {
                RPRM[i] = 0.0;  IPRM[i] = 0;
            }
            int METHOD = 1;          int NDV = 2;           int NCON = 2;
            for (int i = 0; i < NDV; i++)
            {
                X[i] = 1.0; XL[i] = 0.1;    XU[i] = 100.0;
            }
            int IPRINT = 3;            int MINMAX = -1;            int INFO = 0;
            while (true)
            {
                //FortranMethod.DOT(ref INFO, ref METHOD, ref IPRINT, ref NDV, ref NCON, ref X, ref XL, ref XU, ref OBJ, ref MINMAX, ref G, ref RPRM, ref IPRM, ref WK, ref NRWK, ref IWK, ref NRIWK);

                //FortranMethod.DOT(INFO, METHOD, IPRINT, NDV, NCON, X, XL, XU, OBJ, MINMAX, G, RPRM, IPRM, WK, NRWK, IWK, NRIWK);

                FortranMethod.DOT(ref INFO, ref METHOD, ref IPRINT, ref NDV, ref NCON, ref X[0], ref XL[0], ref XU[0], ref OBJ, ref MINMAX, ref G[0], ref RPRM[0], ref IPRM[0], ref WK[0], ref NRWK, ref IWK[0], ref NRIWK);
                //FortranMethod.DOT( NCON,  MINMAX,  NRWK,  NRIWK, ref OBJ, ref X[0], ref XL[0], ref XU[0], ref G[0], ref RPRM[0], ref IPRM[0], ref WK[0], ref IWK[0]);
                //FortranMethod.TOT(ref MINMAX, ref NCON, ref NRWK, ref NDV, ref INFO, ref NRIWK, ref OBJ, ref X[0], ref XL[0], ref XU[0], ref G[0], ref RPRM[0], ref IPRM[0], ref WK[0], ref IWK[0]);
                //FortranMethod.DOT(ref INFO, ref METHOD, ref IPRINT, ref NDV, ref NCON, MINMAX, NRWK, NRIWK, ref OBJ, ref X[0], ref XL[0], ref XU[0], ref G[0], ref RPRM[0], ref IPRM[0], ref WK[0], ref IWK[0]);
                if (INFO == 0)
                {
                    break;
                }
                FortranMethod.EVAL(ref OBJ, ref X[0], ref G[0]);
            }

        }
    }
    public static class FortranMethod
    {
        [DllImport("DOTBOX.dll", CallingConvention = CallingConvention.StdCall)]
        //[DllImport("RPRM.dll", CallingConvention = CallingConvention.StdCall)]
        public static extern double DOT(ref int INFO, ref int METHOD, ref int IPRINT, ref int NDV, ref int NCON, ref double X, ref double XL, ref double XU, ref double OBJ,ref int MINMAX, ref double G, ref double RPRM, ref int IPRM, ref double WK, ref int NRWK, ref int IWK, ref int NRIWK);
        [DllImport("EVAL.dll", CallingConvention = CallingConvention.StdCall)]
        public static extern double EVAL(ref double A, ref double B, ref double C);
        //public static extern double DOT( int NCON,  int MINMAX,  int NRWK,  int NRIWK, ref double OBJ, ref double X, ref double XL, ref double XU, ref double G, ref double RPRM, ref int IPRM, ref double WK, ref int IWK);
        //public static extern double TOT( ref int MINMAX, ref int NCON, ref int NRWK, ref int NDV, ref int INFO, ref int NRIWK, ref double OBJ, ref double X, ref double XL, ref double XU, ref double G, ref double RPRM, ref int IPRM, ref double WK, ref int IWK);
        //public static extern double DOT(ref int INFO, ref int METHOD, ref int IPRINT, ref int NDV, ref int NCON, int MINMAX, int NRWK, int NRIWK, ref double OBJ, ref double X, ref double XL, ref double XU, ref double G, ref double RPRM, ref int IPRM, ref double WK, ref int IWK);
    }



    //[DllImport("DOOT1.dll")]
    //public static extern double DOT(ref int INFO, ref int METHOD, ref int IPRINT, ref int NDV, ref int NCON, ref double[] X, ref double[] XL, ref double[] XU, ref double OBJ, ref int MINMAX, ref double[] G, ref double[] RPRM, ref int[] IPRM, ref double[] WK, ref int NRWK, ref int[] IWK, ref int NRIWK);
    //public static extern double DOT( int INFO,  int METHOD,  int IPRINT,  int NDV,  int NCON, double[] X,  double[] XL, double[] XU,  double OBJ,  int MINMAX,  double[] G, double[] RPRM, int[] IPRM,  double[] WK, int NRWK,  int[] IWK, int NRIWK);
    //

}

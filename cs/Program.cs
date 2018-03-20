using System;
using System.Diagnostics;
using Util;

namespace polyRoots
{
    public class Program
    {
        private static void test()
        {
            Poly2 poly = new Poly2(
                    new Complex(2, 3),
                    new Complex(1, -2),
                    new Complex(-1, 4));
            poly.solve();
            poly.print();
        }

        private static void test1()
        {
            Poly3 poly = new Poly3(
                    new Complex(2, 3),
                    new Complex(1, -2),
                    new Complex(-1, 4),
                    new Complex(5, 6));
            poly.solve();
            poly.print();
        }

        private static void test2()
        {//x^4 - 10·x^3 + 35·x^2 - 50·x + 24 has 1,2,3,4
            double a = 1, b = -10, c = 35, d = -50, e = 24;
            Poly4 poly = new Poly4(
                    new Complex(a, 0),
                    new Complex(b, 0),
                    new Complex(c, 0),
                    new Complex(d, 0),
                    new Complex(e, 0));
            poly.solve();
            poly.print();
        }

        private static void test2a()
        {
            Poly4 poly = new Poly4(
                    new Complex(1.5, -0.5),
                    new Complex(2.1, 0.5),
                    new Complex(-1.2, -0.4),
                    new Complex(-1.6, 0.4),
                    new Complex(2.5, 0.2));
            /*x = -1.5031210042193886368 - 1.1996324925961071257· ∨
             * x = 0.75088813963150753935 - 0.49378327420458781439· ∨ x = 0.62146832688363632989 + 0.53592730551193332322·
             * x = -1.0292354622962654483 + 0.43748846128867246989·*/
            poly.solve
                ();
            poly.check();
            poly.print();
        }

        private static void test2b()
        {
            Poly4 poly = new Poly4(
                    new Complex(1.5, 0),
                    new Complex(2.1, 0),
                    new Complex(-1.2, 0),
                    new Complex(-1.6, 0),
                    new Complex(2.5, 0));
            poly.solve();
            poly.check();
            poly.print();
            /*
            x = 0.68570649639504676818 - 0.51544340371595721844· ∨ x = 0.68570649639504676818 + 0.51544340371595721844· ∨ x = -1.3857064963948774975 - 0.58710645620335704451· ∨ x = -1.3857064963948774975 + 0.58710645620335704451·
             */
        }

        private static void testBook()
        {
            Poly4 poly = new Poly4(
                    new Complex(1, 0),
                    new Complex(6, 0),
                    new Complex(-8, 0),
                    new Complex(-22, 0),
                    new Complex(-105, 0));
            poly.solve();
            poly.check();
            poly.print();
        }


        private static void testBookError()
        {
            Poly4 poly = new Poly4(
                    new Complex(-2, -2),
                    new Complex(1.5, -0.5),
                    new Complex(0, 0.25),
                    new Complex(0, 0),
                    new Complex(0.75, -1.75));
            poly.solve();
            poly.check();
            poly.print();
        }

        static MT19937 rand = new MT19937(0);
        private static double r()
        {
            int r = (int)(rand.Next() % 16);
            double ret = (r - 8) * 0.25;
            //Console.WriteLine(ret);
            return ret;
        }

        Complex c()
        {
            double re = r();
            double im = r();
            return Complex.createComplex(re, im);
        }

        private static void testRandom()
        {
            Stopwatch sw = new Stopwatch();
            sw.Start();
            for (int i = 0; i < 1000000; i++)
            {
                Complex a = new Complex(r(), r());
                if (a.Equals(Complex.ZERO)) a = Complex.ONE;
                Poly4 poly = new Poly4(
                        a,
                        new Complex(r(), r()),
                        new Complex(r(), r()),
                        new Complex(r(), r()),
                        new Complex(r(), r()));
                //poly.print();
                try
                {
                    poly.solve();
                    poly.check();
                }
                catch {
                    Console.WriteLine("A={0} B={1} C={2} D={3} E={4}", poly.A, poly.B, poly.C, poly.D, poly.E);
                }
            }
            sw.Stop();
            Console.WriteLine("Time={0}", sw.Elapsed);
        }

        public static void Main()
        {
            //badRoot();
            testRandom();
            //testPhp();
            //testJava();
            //testBookError();
        }
    }
}

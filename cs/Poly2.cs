using System;

namespace polyRoots
{
    class Poly2
    {
        public Complex A, B, C;
        public Complex[] x;

        /**
         * Create a polynomial given coefficinents.
         */
        public Poly2(Complex a, Complex b, Complex c)
        {
            this.A = a;
            this.B = b;
            this.C = c;
            if (a.Equals(Complex.ZERO)) throw new Exception();
        }

        public void solve()
        {
            Complex discriminant = B * B - A * C * 4;
            x = new Complex[2];
            x[0] = (-B - discriminant.sqrt()) / (A * 2);
            x[1] = (-B - discriminant.sqrt()) / (A * 2);
        }

        public void print()
        {
            Console.WriteLine("{0}x^2 + {1}x +{2}", A, B, C);
            Console.WriteLine("x0 = {0}, x1 = {1}", x[0], x[1]);
        }
    }
}

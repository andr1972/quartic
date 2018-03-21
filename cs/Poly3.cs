using System;
using Util;

namespace polyRoots
{
    class Poly3
    {
        public Complex A, B, C, D;
        public Complex[] x;
        protected Complex p, q;

        /**
         * Create a polynomial from given coefficients.
         *
         */
        public Poly3(Complex a, Complex b, Complex c, Complex d)
        {
            this.A = a;
            this.B = b;
            this.C = c;
            this.D = d;
            if (a.Equals(Complex.ZERO)) throw new Exception();
        }

        private void canonicalForm()
        {
            p = C / A - B * B / (A * A * 3);
            q = (B * B * B * 2) / (A * A * A * 27) + D / A - B * C / (A * A * 3);
        }

        //x = y - b/3a
        private Complex ytox(Complex y)
        {
            return y - B / (A * 3);
        }

        int find_m_from_e(Complex z, Complex[] e)
        {
            for (int i = 0; i < 3; i++)
            {
                Complex r = z - e[i];
                if (r.abs() < 1e-7)
                {
                    return i;
                }
            }
            throw new Exception(); //it is possible only due machine errors
        }

        int find_l_for_m(int m)
        {
            switch (m)

            {
                case 0: return 0;
                case 1: return 2;
                case 2: return 1;
                default: throw new Exception(); //it is not possible
            }
        }

        public void solve()
        {
            canonicalForm();
            x = new Complex[3];
            if (p.abs() < 1e-12)
            {
                Complex []roots = Complex.ZERO.subtract(q).nthRoots(3);
                for (int i = 0; i < 3; i++)
                    x[i] = ytox(roots[i]);
            }
            else
            {
                Complex z2 = (q * q + p * p * p * 4.0 / 27.0);
                if (z2.abs() < 1e-13) z2 = Complex.ZERO; //double root
                Complex z = z2.sqrt();
                Complex z0 = (z - q) / 2;
                Complex v0 = z0.nthRoot(3);
                Complex u = (-q - z0).nthRoot(3);
                z = v0 * u / p * (-3);
                Complex[] e = new Complex[3];
                e[0] = new Complex(1, 0);
                e[1] = new Complex(-1.0 / 2, Math.Sqrt(3) / 2);
                e[2] = new Complex(-1.0 / 2, -Math.Sqrt(3) / 2);
                int m = find_m_from_e(z, e);
                int n = find_l_for_m(m);
                Complex u0 = e[n] * u;
                Complex[] y = new Complex[3];
                y[0] = v0 + u0;
                y[1] = e[1] * v0 + e[1] * e[1] * u0;
                y[2] = e[1] * e[1] * v0 + e[1] * u0;
                for (int i = 0; i < 3; i++)
                    x[i] = ytox(y[i]);
            }
        }

        public void check()
        {
            for (int i = 0; i < 3; i++)
            {
                Complex xi2 = x[i] * x[i];
                Complex z = A * xi2 * x[i] + B * xi2 + C * x[i] + D;
                if (z.abs() > 1e-9)
                {
                    throw new Exception("bad root");
                }
            }
        }

        public void print()
        {
            Console.WriteLine("{0}x^3 + {1}x^2 +{2}x + {3}", A, B, C, D);
            for (int i = 0; i < 3; i++)
                Console.WriteLine("x{0}= {1}", i, x[i]);
        }
    }
}

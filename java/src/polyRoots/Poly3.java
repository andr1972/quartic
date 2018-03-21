package polyRoots;

import util.Complex;

public class Poly3 {
	public Complex A,B,C,D;
	public Complex[] x;
	public Complex p,q;

	/**
     * Create a polynomial from given coefficients.
	 * @throws Exception
     *
     */
    public Poly3(Complex a, Complex b, Complex c, Complex d) throws Exception {
        this.A = a;
        this.B = b;
        this.C = c;
        this.D = d;
        if (a.equals(Complex.ZERO)) throw new Exception();
    }

    // p = c/a-b*b/(3*a*a);
	// q = 2*b*b*b/(27*a*a*a)+d/a-b*c/(3*a*a);
	private void canonicalForm()
	{
		p = C.divide(A).subtract(B.multiply(B).divide(A.multiply(A).multiply(3)));
		q = (B.multiply(B).multiply(B).multiply(2).divide(A.multiply(A).multiply(A).multiply(27)))
			.add(D.divide(A)).subtract(B.multiply(C).divide(A.multiply(A).multiply(3)));
	}

	//x = y - b/3a
	private Complex ytox(Complex y)
	{
		return y.subtract(B.divide(A.multiply(3)));
	}

	int find_m_from_eps(Complex z, Complex[] eps) throws Exception
	{
		for (int i=0; i<3; i++)
		{
			Complex r = z.subtract(eps[i]);
			if (r.abs()<1e-7)
			{
				return i;
			}
		}
		throw new Exception(); //it is possible only due machine errors
	}

	int find_l_for_m(int m) throws Exception
	{
		switch (m)
		{
		case 0: return 0;
		case 1: return 2;
		case 2: return 1;
		default: throw new Exception(); //it is not possible
		}
	}

	public void solve() throws Exception
	{
		canonicalForm();
		x = new Complex[3];
		if (p.abs()<1e-12)
		{
			Complex roots[] = Complex.ZERO.subtract(q).nthRoots(3);
			for (int i=0; i<3; i++)
				x[i] = ytox(roots[i]);
		}
		else
		{
			Complex t2 = q.multiply(q).add(p.multiply(p).multiply(p).multiply((double)4/27));
			if (t2.abs() < 1e-13) t2 = Complex.ZERO; //double root
			Complex temp = t2.sqrt();
			Complex z0 = Complex.ZERO.subtract(q).add(temp).divide(2);
			Complex v0 = z0.nthRoot(3);
			Complex u = Complex.ZERO.subtract(q).subtract(z0).nthRoot(3);

			temp = v0.multiply(u).divide(p).multiply(-3);
			Complex[] eps = new Complex[3];
			eps[0] = new Complex(1,0);
			eps[1] = new Complex(-1.0/2,Math.sqrt(3)/2);
			eps[2] = new Complex(-1.0/2,-Math.sqrt(3)/2);
			int m = find_m_from_eps(temp,eps);
			int n = find_l_for_m(m);
			Complex u0 = eps[n].multiply(u);
			Complex []y = new Complex[3];
			y[0] = v0.add(u0);
			y[1] = (eps[1].multiply(v0)).add(eps[1].multiply(eps[1]).multiply(u0));
			y[2] = (eps[1].multiply(eps[1]).multiply(v0)).add(eps[1].multiply(u0));
			for (int i=0; i<3; i++)
				x[i] = ytox(y[i]);
		}
	}

	double rootError = 0;

	public void check() throws Exception
    {
		rootError = 0;
        for (int i = 0; i < 3; i++)
        {
            Complex xi2 = x[i].multiply(x[i]);
            Complex z = A.multiply(xi2).multiply(x[i]).
            		add(B.multiply(xi2)).add(C.multiply(x[i])).add(D);
            rootError = Math.max(rootError, z.abs());
            if (z.abs() > 1e-9)
            {
                //Console.WriteLine("{0}x^4 + {1}x^3 + {2}x^2 + {3}x + {4}", A, B, C, D, E);
            	throw new Exception("bad root");
                //break;
            }
        }
    }

	public void print()
    {
    	System.out.format("%sx^3 +%sx^2 + %sx +%s\n",
    			A.toString(),B.toString(),
    			C.toString(),D.toString());
    }
}

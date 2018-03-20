package polyRoots;

import util.Complex;

public class Poly4 {
	public Complex A,B,C,D,E;
	public Complex[] x;
	public Complex p,q,r;

	/**
     * Create a polynomial from given coefficients.
	 * @throws Exception
     *
     */
    public Poly4(Complex a, Complex b, Complex c, Complex d, Complex h) throws Exception {
        this.A = a;
        this.B = b;
        this.C = c;
        this.D = d;
        this.E = h;
        if (a.equals(Complex.ZERO)) throw new Exception();
    }

    Complex find_maxB2_4A(Complex[] cubicRoots)
    {
    	int k = 0;
    	double maxModulus = 0;
		for (int i=0; i<3; i++)
		{
			Complex y = cubicRoots[i];
			Complex v = B.multiply(B).subtract(A.multiply(y).multiply(4));
			double modulus = v.abs();
			if (modulus>maxModulus)
			{
				maxModulus = modulus;
				k = i;
			}
		}
		return cubicRoots[k];
    }

    //Neumark - Solution Of Cubic & Quartic Equations
    public void solve() throws Exception
    {
        x = new Complex[4];
        if (B.abs() < 1e-12 && D.abs() < 1e-12)
        {
            Poly2 poly2 = new Poly2(A, C, E);
            poly2.solve();
            x[0] = poly2.x[0].sqrt();
            x[1] = x[0].negate();
            x[2] = poly2.x[1].sqrt();
            x[3] = x[2].negate();
        }
        else
        {
        	Poly3 subpoly = new Poly3(Complex.ONE,
            		C.multiply(-2),
            		C.multiply(C).add(B.multiply(D)).subtract(A.multiply(E).multiply(4)),
            		B.multiply(C).multiply(D).subtract(B.multiply(B).multiply(E)).subtract(A.multiply(D).multiply(D)).negate());
            subpoly.solve();
            Complex y = find_maxB2_4A(subpoly.x);
            Complex B2_4A = B.multiply(B).subtract(A.multiply(y).multiply(4)).sqrt();
            Complex G,g,H,h;
            if (B2_4A.abs()<1e-12)
            {
            	G = B.multiply(0.5);
	            g = G;
	            H = C.subtract(y).multiply(0.5);
	            h = H;
            }
            else
            {
	            G = B.add(B2_4A).multiply(0.5);
	            g = B.subtract(B2_4A).multiply(0.5);
	            Complex part = B.multiply(C.subtract(y)).subtract(A.multiply(D).multiply(2)).divide(B2_4A).multiply(0.5);
	            H = C.subtract(y).multiply(0.5).add(part);
	            h = C.subtract(y).multiply(0.5).subtract(part);
            }
            Poly2 poly2a = new Poly2(A, G, H);
            poly2a.solve();
            Poly2 poly2b = new Poly2(A, g, h);
            poly2b.solve();
            x[0] = poly2a.x[0];
            x[1] = poly2a.x[1];
            x[2] = poly2b.x[0];
            x[3] = poly2b.x[1];
        }
    }

    public double rootError = 0;

    public void check() throws Exception
    {
    	rootError = 0;
        for (int i = 0; i < 4; i++)
        {
            Complex xi2 = x[i].multiply(x[i]);
            Complex z = A.multiply(xi2).multiply(xi2).add(B.multiply(xi2).multiply(x[i])).
            		add(C.multiply(xi2)).add(D.multiply(x[i])).add(E);
            rootError = Math.max(rootError, z.abs());
            if (z.abs() > 1e-9)
            {
            	for (int j = 0; j < 4; j++)
            	{
            		System.out.println("x"+j+ ":" + x[j]);
            	}
            	System.out.println(A.toString() + "x^4+" + B.toString()+"x^3 + "+C.toString()+
            			"x^2"+D.toString()+"x + "+E.toString());
            	throw new Exception("bad root");
            }
        }
    }


	public void print()
    {
		System.out.format("%sx^4 + %sx^3 +%sx^2 + %sx +%s\n",
    			A.toString(),B.toString(),
    			C.toString(),D.toString(),
    			E.toString());
    }
}

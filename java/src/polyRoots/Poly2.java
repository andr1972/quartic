package polyRoots;

import util.Complex;

public class Poly2 {
	public Complex A,B,C;
	public Complex[] x;

	/**
     * Create a polynomial given coefficinents.
	 * @throws Exception
     *
     */
    public Poly2(Complex a, Complex b, Complex c) throws Exception {
        this.A = a;
        this.B = b;
        this.C = c;
        if (a.equals(Complex.ZERO)) throw new Exception();
    }

    public void solve()
    {
    	Complex discriminant = B.multiply(B).subtract(A.multiply(C).multiply(4));
    	x = new Complex[2];
    	x[0] = Complex.ZERO.subtract(B).subtract(discriminant.sqrt()).divide(A.multiply(2));
    	x[1] = Complex.ZERO.subtract(B).add(discriminant.sqrt()).divide(A.multiply(2));
    }

    public void print()
    {
    	System.out.format("(%f,%f)x^2 + (%f,%f)x +(%f,%f)\n",
    			A.getReal(),A.getImaginary(),
    			B.getReal(),B.getImaginary(),
    			C.getReal(),C.getImaginary());
    	System.out.format("x0 = (%f,%f), x1 = (%f,%f)\n",
    			x[0].getReal(),x[0].getImaginary(),
    			x[1].getReal(),x[1].getImaginary());
    }
}

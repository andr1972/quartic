package polyRoots;
import java.util.Locale;

import util.Complex;
import util.MersenneTwisterFast;

//trzeba arbitrary precision
public class Main {

	private static void test1()
	{
	    try {
			Poly3 poly = new Poly3(
					new Complex(2,3),
					new Complex(1,-2),
					new Complex(-1,4),
					new Complex(5,6));
			poly.solve();
			poly.print();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private static void test1a()
	{
	    try {
			Poly3 poly = new Poly3(
					new Complex(1,0),
					new Complex(5,0),
					new Complex(4,0),
					new Complex(0,0));
			poly.solve();
			poly.print();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private static void test2()
	{//x^4 - 10·x^3 + 35·x^2 - 50·x + 24 has 1,2,3,4
	    try {
	    	double a=1, b=10, c=35, d=-50, h=24;
			Poly4 poly = new Poly4(
					new Complex(a,0),
					new Complex(b,0),
					new Complex(c,0),
					new Complex(d,0),
					new Complex(h,0));
			poly.solve();
			poly.print();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private static void test2a()
	{ //bi square x^4 - 3·x^2 + 2
	    try {
	    	double a=1, b=0, c=-3, d=0, e=2;
			Poly4 poly = new Poly4(
					new Complex(a,0),
					new Complex(b,0),
					new Complex(c,0),
					new Complex(d,0),
					new Complex(e,0));
			poly.solve();
			poly.print();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private static void test()
	{
	    try {
			Poly2 poly = new Poly2(
					new Complex(2,3),
					new Complex(1,-2),
					new Complex(-1,4));
			poly.solve();
			poly.print();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	static MersenneTwisterFast rand = new MersenneTwisterFast(0);

	private static double r() {
		double ret = (rand.next() % 16 - 8)*0.25;
		return ret;
	}

	private static Complex c()
	{
		double re = r();
		double im = r();
		return Complex.createComplex(re,im);
	}


	private static void test3random() {
		double max = 0;
		for (int i=0; i<3; i++)
		{
			try {
				Complex a = new Complex(r(),r());
				if (a.equals(Complex.ZERO))
					a = Complex.ONE;
				Poly3 poly = new Poly3(
						a,
						new Complex(r(),r()),
						new Complex(r(),r()),
						new Complex(r(),r()));
				poly.print();
				poly.solve();
				poly.check();
				max = Math.max(max, poly.rootError);

			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		System.out.println("rootError=" + max);
	}

	private static void test4random() {
		double max = 0;
		long startTime = System.nanoTime();
		for (int i=0; i<1000000; i++)
		{
			try {
				Complex A = c();
				if (A.equals(Complex.ZERO))
					A = Complex.ONE;
				Complex B = c();
				Complex C = c();
				Complex D = c();
				Complex E = c();
				Poly4 poly = new Poly4(A,B,C,D,E);
				poly.solve();
				//poly.check();
				//poly.print();
				max = Math.max(max, poly.rootError);

			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		long estimatedTime = System.nanoTime() - startTime;
		System.out.println("Time=" + estimatedTime/1e9);
	}

	//test 3 roots
	private static void roots3()
	{ //bi square x^4 - 3·x^2 + 2
	    try {
	    	Complex a = Complex.ONE;
	    	Complex b = Complex.ONE;
	    	Complex c = new Complex(- 3.0/4.0, -1.0/8.0);
	    	Complex d = new Complex(Math.sqrt(6*Math.sqrt(130) + 66)/16.0 - 1/2.0, Math.sqrt(6*Math.sqrt(130) - 66)/16.0 - 1/16.0);
	    	Complex e = new Complex(Math.sqrt(6*Math.sqrt(130) + 66)/64.0 - 53/256.0, Math.sqrt(6*Math.sqrt(130) - 66)/64.0 - 1/32.0);

			Poly4 poly = new Poly4(
					a,
					b,
					c,
					d,
					e);
			poly.solve();
			poly.check();
			poly.print();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/*
	 * x0=  1.1 - 1.1·i
	   x1 = -1.1 + 0.2·i
	 * */
	private static void doubleRoots()
	{ //bi square x^4 - 3·x^2 + 2
	    try {
	    	Complex a = Complex.ONE;
	    	Complex b = new Complex(0, 1.8);
	    	Complex c = new Complex(-2.79, 2.86);
	    	Complex d = new Complex(-2.574, -1.782);
	    	Complex e = new Complex(-1.0648, -2.8314);

			Poly4 poly = new Poly4(
					a,
					b,
					c,
					d,
					e);
			poly.solve();
			poly.check();
			poly.print();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private static void testBad()
	{
	    try {
	    	Complex a = new Complex(-0.25, 0.0);
	    	Complex b = new Complex(0.0, -0.5);
	    	Complex c = new Complex(0.75, -0.75);
	    	Complex d = new Complex(0.75, 0.5);
	    	Complex e = new Complex(0.5, 0.5);

			Poly4 poly = new Poly4(
					a,
					b,
					c,
					d,
					e);
			poly.solve();
			poly.check();
			poly.print();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private static void quadruple()
	{
	    try {
	    	Complex a = new Complex(1,0);
	    	Complex b = new Complex(2,0);
	    	Complex c = new Complex(1.5,0);
	    	Complex d = new Complex(0.5,0);
	    	Complex e = new Complex(1/16.0,0);

			Poly4 poly = new Poly4(
					a,
					b,
					c,
					d,
					e);
			poly.solve();
			poly.check();
			poly.print();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private static void quadrupleComplex()
	{
	    try {
	    	Complex a = new Complex(1,0);
	    	Complex b = new Complex(-4,-4);
	    	Complex c = new Complex(0,12);
	    	Complex d = new Complex(8,-8);
	    	Complex e = new Complex(-4.0);

			Poly4 poly = new Poly4(
					a,
					b,
					c,
					d,
					e);
			poly.solve();
			poly.check();
			poly.print();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private static void badE()
	{
	    try {
			Complex a = new Complex(-0.750000, 1.250000);
            Complex b = new Complex(-2.000000, 1.250000);
            Complex c = new Complex(0.500000, 0.500000);
            Complex d = new Complex(-0.250000, -2.000000);
            Complex e = new Complex(0.250000, -0.750000);

			Poly4 poly = new Poly4(
					a,
					b,
					c,
					d,
					e);
			poly.solve();
			//poly.check();
			//poly.print();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		Locale.setDefault(new Locale("en", "US"));
		//roots3();
		//doubleRoots();
		//testBad();
		test4random();
		//badE();
		//quadruple();
		//quadrupleComplex();
	}
}

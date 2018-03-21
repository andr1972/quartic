#include "Poly4.h"
#include "Poly2.h"
#include "Poly3.h"

/**
* Create a polynomial from given coefficients.
*
*/
Poly4::Poly4(Complex &a, Complex &b, Complex &c, Complex &d, Complex &h)
{
	this->A = a;
	this->B = b;
	this->C = c;
	this->D = d;
	this->E = h;
	if (a.Equals(Complex::ZERO)) throw "A can't be zero";
}

Complex Poly4::findMaxDelta(Complex* cubicRoots)
{
	int k = 0;
	double maxModulus = 0;
	for (int i = 0; i < 3; i++)
	{
		Complex y = cubicRoots[i];
		Complex Delta = B*B - A*y*4;
		double modulus = Delta.abs();
		if (modulus > maxModulus)
		{
			maxModulus = modulus;
			k = i;
		}
	}
	return cubicRoots[k];
}

//Neumark - Solution Of Cubic & Quartic Equations
void Poly4::solve()
{
	if (B.abs() < 1e-12 && D.abs() < 1e-12)
	{
		Poly2 poly2(A, C, E);
		poly2.solve();
		x[0] = poly2.x[0].sqrt();
		x[1] = -x[0];
		x[2] = poly2.x[1].sqrt();
		x[3] = -x[2];
	}
	else
	{
		Poly3 subpoly(Complex::ONE,
			C*(-2),
			C*C + B*D - A*E*4,
			-(B*C*D - B*B*E - A*D*D));
		subpoly.solve();
		Complex y = findMaxDelta(subpoly.x);
		Complex DeltaRoot = (B*B - A*y*4).sqrt();
		Complex G, g, H, h;
		if (DeltaRoot.abs() < 1e-12)
		{
			G = B*0.5;
			g = G;
			H = (C - y)* 0.5;
			h = H;
		}
		else
		{
			G = (B + DeltaRoot)* 0.5;
			g = (B - DeltaRoot)* 0.5;
			Complex part = (B* (C-y) - A*D*2)/(DeltaRoot*2);
			H = (C - y)* 0.5 + part;
			h = (C - y)* 0.5 - part;
		}
		Poly2 poly2a(A, G, H);
		poly2a.solve();
		Poly2 poly2b(A, g, h);
		poly2b.solve();
		x[0] = poly2a.x[0];
		x[1] = poly2a.x[1];
		x[2] = poly2b.x[0];
		x[3] = poly2b.x[1];
	}
}

void Poly4::check()
{
	for (int i = 0; i < 4; i++)
	{
		Complex xi2 = x[i]* x[i];
		Complex z = A*xi2*xi2 + B*xi2*x[i] + C*xi2 + D*x[i] + E;
		if (z.abs() > 1e-9)
		{
			throw "bad root";
		}
	}
}


void Poly4::print()
{
	printf("%sx^4 + %sx^3 + %sx^2 + %sx + %s\n", A.toString().c_str(), B.toString().c_str(),
		C.toString().c_str(), D.toString().c_str(), E.toString().c_str());
}

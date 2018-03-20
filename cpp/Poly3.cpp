#include "Poly3.h"

/**
* Create a polynomial from given coefficients.
*
*/
Poly3::Poly3(Complex &a, Complex &b, Complex &c, Complex &d)
{
	this->A = a;
	this->B = b;
	this->C = c;
	this->D = d;
	if (a.Equals(Complex::ZERO)) throw "A can't be 0";
}

void Poly3::canonicalForm()
{
	p = C/A - B*B/(A*A*3);
	q = (B*B*B*2)/(A*A*A*27) + D/A - B*C/(A*A*3);
}

Complex Poly3::ytox(Complex y)
{
	return y - B/(A*3);
}

int Poly3::find_m_from_eps(Complex z, Complex *eps)
{
	for (int i = 0; i < 3; i++)
	{
		Complex r = z - eps[i];
		if (r.abs() < 1e-7)
		{
			return i;
		}
	}
	throw "not found e"; //it is possible only due machine errors
}

int Poly3::find_l_for_m(int m)
{
	switch (m)
	{
	case 0: return 0;
	case 1: return 2;
	case 2: return 1;
	default: throw "is not possible";
	}
}

void Poly3::solve()
{
	canonicalForm();
	if (p.abs() < 1e-12)
	{
		vector<Complex> roots = (-q).nthRoots(3);
		for (int i = 0; i < 3; i++)
			x[i] = ytox(roots[i]);
	}
	else
	{
		Complex temp = (q* q + p*p*p* 4.0 / 27.0).sqrt();
		Complex z0 = temp-q;
		z0 = z0 / 2;
		Complex v0 = z0.nthRoot(3);
		Complex u = (-q-z0).nthRoot(3);
		temp = v0 * u / p* (-3);
		Complex eps[3];
		eps[0] = Complex::createComplex(1, 0);
		eps[1] = Complex::createComplex(-1.0 / 2, sqrt(3) / 2);
		eps[2] = Complex::createComplex(-1.0 / 2, -sqrt(3) / 2);
		int m = find_m_from_eps(temp, eps);
		int n = find_l_for_m(m);
		Complex u0 = eps[n]* u;
		Complex y[3];
		y[0] = v0 + u0;
		y[1] = eps[1]* v0 + eps[1]* eps[1]* u0;
		y[2] = eps[1]* eps[1]* v0 + eps[1]* u0;
		for (int i = 0; i < 3; i++)
			x[i] = ytox(y[i]);
	}
}

void Poly3::check()
{
	for (int i = 0; i < 3; i++)
	{
		Complex xi2 = x[i]* x[i];
		Complex z = A* xi2* x[i] + B * xi2 + C * x[i] + D;
		if (z.abs() > 1e-9)
		{
			throw "bad root";
		}
	}
}

void Poly3::print()
{
}
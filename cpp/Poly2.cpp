#include "Poly2.h"

Poly2::Poly2(Complex &a, Complex &b, Complex &c)
{
	this->A = a;
	this->B = b;
	this->C = c;
	if (a.Equals(Complex::ZERO)) throw "a must be != 0";
}

Poly2::~Poly2()
{
}

void Poly2::solve()
{
	Complex discriminant = B*B - A*C*4;
	x[0] = (-B - discriminant.sqrt()) / (A*2);
	x[1] = (-B + discriminant.sqrt()) / (A*2);
}

void Poly2::print()
{
	printf("%sx^2 + %sx + %s\n", A.toString().c_str(), B.toString().c_str(), C.toString().c_str());
	printf("x0 = %s, x1 = %s\n", x[0].toString().c_str(), x[1].toString().c_str());
}
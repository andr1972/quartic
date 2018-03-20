#pragma once
#include "Complex.h"

class Poly2
{
public:
	Complex A, B, C;
	Complex x[2];
	Poly2(Complex &a, Complex &b, Complex &c);
	~Poly2();
	void solve();
	void print();
};


#pragma once
#include "Complex.h"

class Poly3
{
	Complex p, q;
	void canonicalForm();
	Complex ytox(Complex y);
	int find_m_from_eps(Complex z, Complex *eps);
	int find_l_for_m(int m);
public:
	Complex A, B, C, D;
	Complex x[3];
	Poly3(Complex &a, Complex &b, Complex &c, Complex &d);
	~Poly3(){};
	void solve();
	void check();
	void print();
};


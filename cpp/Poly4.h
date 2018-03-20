#pragma once
#include "Complex.h"

class Poly4
{
	Complex find_maxB2_4A(Complex* cubicRoots);
public:
	Complex A, B, C, D, E;
	Complex x[4];
	Poly4(Complex &a, Complex &b, Complex &c, Complex &d, Complex &h);
	~Poly4(){};
	void solve();
	void check();
	void print();
};


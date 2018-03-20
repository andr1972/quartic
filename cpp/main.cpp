#include "Complex.h"
#include "Poly4.h"
#include "Poly3.h"
#include <random>
#include <chrono>

using namespace std;


mt19937 gen(0);

double r()
{
	int r = (int)(gen() % 16);
	double ret = (r - 8) * 0.25;
	//printf("%f\n",ret);
	return ret;
}

Complex c()
{
	double re = r();
	double im = r();
	return Complex::createComplex(re,im);
}

void testBad()
{
	Poly4 poly(
		Complex::createComplex(-1.250000, -0.250000),
		Complex::createComplex(-0.750000, -0.500000),
		Complex::createComplex(0.500000, 2.000000),
		Complex::createComplex(-2.000000, -1.500000),
		Complex::createComplex(2.000000, 1.000000)	);
	poly.solve();
	poly.check();
	poly.print();
}

void testRandom()
{
	auto time_a = chrono::high_resolution_clock::now();
	for (int i = 0; i < 1000000; i++)
	{
		Complex A = c();
		if (A.Equals(Complex::ZERO)) A = Complex::ONE;
		Complex B = c();
		Complex C = c();
		Complex D = c();
		Complex E = c();
		Poly4 poly(A,B,C,D,E);
		//poly.print();
		try
		{
			poly.solve();
			poly.check();
		}
		catch (...){
			printf("A=%s B=%s C=%s D=%s E=%s\n", poly.A.toString().c_str(), poly.B.toString().c_str(),
				poly.C.toString().c_str(), poly.D.toString().c_str(), poly.E.toString().c_str());
		}
	}
	auto time_b = chrono::high_resolution_clock::now();
	chrono::duration<double> elapsed_seconds = time_b - time_a;
	printf("time=%f sec.", elapsed_seconds.count());
}

int main()
{
	testRandom();
	//testBad();
    return 0;
}


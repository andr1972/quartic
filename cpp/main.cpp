#include "Complex.h"
#include "Poly4.h"
#include "Poly3.h"
#include "Poly2.h"
#include <random>
#include <chrono>
#include <algorithm>

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

void test()
{
	Poly3 poly(
		Complex::createComplex(1, 0),
		Complex::createComplex(-1.5, -0.6),
		Complex::createComplex(0.45 , 0.88),
		Complex::createComplex(-0.087, -0.57));
	poly.solve();
	poly.print();
}

void randomRoots()
{
	vector<double>dist;
	for (int i = 0; i < 100000; i++)
	{
		Complex roots[3];
		double a = r();
		double b = r();
		double c = r();
		double d = r();
		double e = r();
		double f = r();
		roots[0] = Complex::createComplex(a,b);
		roots[1] = Complex::createComplex(c,d);
		roots[2] = Complex::createComplex(e,f);
		Complex A(1,0);
		Complex B(-a-c-e, -b-d-f);
		Complex C(a*c+a*e+c*e-b*d-b*f-d*f, b*c+a*d+b*e+d*e+a*f+c*f);
		Complex D(-a*c*e + b*d*e+ b*c*f+a*d*f, - b* c* e - a* d* e -  a* c* f + b* d* f);
		Poly3 poly(A, B, C, D);
		poly.solve();
		double maxDist = poly.maxDist(roots);
		if (maxDist / 4e-16>100000000)
			printf("%f\n", maxDist/4e-16);
		dist.push_back(maxDist);
	}
	sort(dist.begin(), dist.end());
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
			//poly.check();
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
	//test();
    return 0;
}


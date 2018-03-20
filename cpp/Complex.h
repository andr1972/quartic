#pragma once
#include <vector>

using namespace std;

class Complex
{
	/** The real part. */
	double real;
	/** The imaginary part. */
	double imaginary;
	/** Record whether this Complex1 number is equal to NaN. */
	bool IsNaN;
	/** Record whether this Complex1 number is infinite. */
	bool IsInfinite;
public:
	static Complex createComplex(double realPart, double imaginaryPart);
	/** The square root of -1. A number representing "0.0 + 1.0i" */
	static Complex I;
	// CHECKSTYLE: stop ConstantName
	/** A Complex1 number representing "NaN + NaNi" */
	static Complex NaN;
	// CHECKSTYLE: resume ConstantName
	/** A Complex1 number representing "+INF + INFi" */
	static Complex INF;
	/** A Complex1 number representing "1.0 + 0.0i" */
	static Complex ONE;
	/** A Complex1 number representing "0.0 + 0.0i" */
	static Complex ZERO;
	Complex();
	Complex(double real);
	Complex(double real, double imaginary);
	~Complex() {};
	Complex operator=(Complex & src);
	double getImaginary(){ return imaginary; };
	double getReal(){ return real; };
	double abs();
	Complex operator+(Complex &addend);
	//Complex add(Complex &addend);
	//Complex add(double addend);
	Complex operator+(double addend);
	Complex conjugate();
	Complex operator/(Complex &divisor);
	Complex operator/(double divisor);
	//Complex divide(Complex &divisor);
	//Complex divide(double divisor);
	Complex reciprocal();
	bool Equals(Complex &c);
	static bool equals(Complex &x, Complex &y, double eps);
	Complex operator*(Complex &factor);
	//Complex multiply(Complex &factor);
	//Complex multiply(int factor);
	Complex operator*(int factor);
	//Complex multiply(double factor);
	Complex operator*(double factor);
	//Complex negate();
	Complex operator-();
	//Complex subtract(Complex &subtrahend);
	Complex operator-(Complex &subtrahend);
	Complex operator-(double subtrahend);
	//Complex subtract(double subtrahend);
	Complex log();
	Complex acos();
	Complex asin();
	Complex atan();
	Complex cos();
	Complex cosh();
	Complex exp();
	Complex pow(Complex &x);
	Complex pow(double x);
	Complex sin();
	Complex sinh();
	Complex sqrt();
	Complex sqrt1z();
	Complex tan();
	Complex tanh();
	double getArgument();
	Complex nthRoot(int n);
	vector<Complex> nthRoots(int n);
	string toString();
};


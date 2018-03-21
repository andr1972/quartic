#include "Complex.h"
#include <math.h>
#include <string>

#define M_PI		3.14159265358979323846

Complex Complex::I(0.0, 1.0);
Complex Complex::NaN(NAN, NAN);
Complex Complex::INF(INFINITY, INFINITY);
Complex Complex::ONE(1.0, 0.0);
Complex Complex::ZERO(0.0, 0.0);

Complex::Complex() : Complex(0.0, 0.0)
{ }

/**
* Create a Complex number given only the real part.
*
* @param real Real part.
*/
Complex::Complex(double real) : Complex(real, 0.0)
{ }

/**
* Create a Complex number given the real and imaginary parts.
*
* @param real Real part.
* @param imaginary Imaginary part.
*/
Complex::Complex(double real, double imaginary)
{
	this->real = real;
	this->imaginary = imaginary;

	IsNaN = isnan(real) || isnan(imaginary);
	IsInfinite = !IsNaN &&
		(isinf(real) || isinf(imaginary));
}

Complex Complex::operator=(const Complex &src)
{
	this->real = src.real;
	this->imaginary = src.imaginary;
	this->IsNaN = src.IsNaN;
	this->IsInfinite = src.IsInfinite;
	return *this;
}

/**
* Return the absolute value of this Complex number.
* Returns {@code NaN} if either real or imaginary part is {@code NaN}
* and {@code Double.POSITIVE_INFINITY} if neither part is {@code NaN},
* but at least one part is infinite.
*
* @return the absolute value.
*/
double Complex::abs()
{
	if (IsNaN)
	{
		return NAN;
	}
	if (IsInfinite)
	{
		return INFINITY;
	}
	if (fabs(real) < fabs(imaginary))
	{
		if (imaginary == 0.0)
		{
			return fabs(real);
		}
		double q = real / imaginary;
		return fabs(imaginary) * ::sqrt(1 + q * q);
	}
	else
	{
		if (real == 0.0)
		{
			return fabs(imaginary);
		}
		double q = imaginary / real;
		return fabs(real) * ::sqrt(1 + q * q);
	}
}

/**
* Returns a {@code Complex} whose value is
* {@code (this + addend)}.
* Uses the definitional formula
* <p>
*   {@code (a + bi) + (c + di) = (a+c) + (b+d)i}
* </p>
* If either {@code this} or {@code addend} has a {@code NaN} value in
* either part, {@link #NaN} is returned; otherwise {@code Infinite}
* and {@code NaN} values are returned in the parts of the result
* according to the rules for {@link java.lang.Double} arithmetic.
*
* @param  addend Value to be added to this {@code Complex}.
* @return {@code this + addend}.
*/
Complex Complex::operator+(const Complex &addend)
{
	if (IsNaN || addend.IsNaN)
	{
		return NaN;
	}
	return createComplex(real + addend.getReal(),
		imaginary + addend.getImaginary());
}

/**
* Returns a {@code Complex} whose value is {@code (this + addend)},
* with {@code addend} interpreted as a real number.
*
* @param addend Value to be added to this {@code Complex}.
* @return {@code this + addend}.
* @see #add(Complex)
*/
Complex Complex::operator+(double addend)
//Complex Complex::add(double addend)
{
	if (IsNaN || isnan(addend))
	{
		return NaN;
	}
	return createComplex(real + addend, imaginary);
}

/**
* Returns the conjugate of this Complex number.
* The conjugate of {@code a + bi} is {@code a - bi}.
* <p>
* {@link #NaN} is returned if either the real or imaginary
* part of this Complex number equals {@code Double.NaN}.
* </p><p>
* If the imaginary part is infinite, and the real part is not
* {@code NaN}, the returned value has infinite imaginary part
* of the opposite sign, e.g. the conjugate of
* {@code 1 + POSITIVE_INFINITY i} is {@code 1 - NEGATIVE_INFINITY i}.
* </p>
* @return the conjugate of this Complex object.
*/
Complex Complex::conjugate()
{
	if (IsNaN)
	{
		return NaN;
	}
	return createComplex(real, -imaginary);
}

/**
* Returns a {@code Complex} whose value is
* {@code (this / divisor)}.
* Implements the definitional formula
* <pre>
*  <code>
*    a + bi          ac + bd + (bc - ad)i
*    ----------- = -------------------------
*    c + di         c<sup>2</sup> + d<sup>2</sup>
*  </code>
* </pre>
* but uses
* <a href="http://doi.acm.org/10.1145/1039813.1039814">
* prescaling of operands</a> to limit the effects of overflows and
* underflows in the computation.
* <p>
* {@code Infinite} and {@code NaN} values are handled according to the
* following rules, applied in the order presented:
* <ul>
*  <li>If either {@code this} or {@code divisor} has a {@code NaN} value
*   in either part, {@link #NaN} is returned.
*  </li>
*  <li>If {@code divisor} equals {@link #ZERO}, {@link #NaN} is returned.
*  </li>
*  <li>If {@code this} and {@code divisor} are both infinite,
*   {@link #NaN} is returned.
*  </li>
*  <li>If {@code this} is finite (i.e., has no {@code Infinite} or
*   {@code NaN} parts) and {@code divisor} is infinite (one or both parts
*   infinite), {@link #ZERO} is returned.
*  </li>
*  <li>If {@code this} is infinite and {@code divisor} is finite,
*   {@code NaN} values are returned in the parts of the result if the
*   {@link java.lang.Double} rules applied to the definitional formula
*   force {@code NaN} results.
*  </li>
* </ul>
*
* @param divisor Value by which this {@code Complex} is to be divided.
* @return {@code this / divisor}.
*/
Complex Complex::operator/(const Complex &divisor)
{
	if (IsNaN || divisor.IsNaN)
	{
		return NaN;
	}

	double c = divisor.getReal();
	double d = divisor.getImaginary();
	if (c == 0.0 && d == 0.0)
	{
		return NaN;
	}

	if (divisor.IsInfinite && !IsInfinite)
	{
		return ZERO;
	}

	if (fabs(c) < fabs(d))
	{
		double q = c / d;
		double denominator = c * q + d;
		return createComplex((real * q + imaginary) / denominator,
			(imaginary * q - real) / denominator);
	}
	else
	{
		double q = d / c;
		double denominator = d * q + c;
		return createComplex((imaginary * q + real) / denominator,
			(imaginary - real * q) / denominator);
	}
}

/**
* Returns a {@code Complex} whose value is {@code (this / divisor)},
* with {@code divisor} interpreted as a real number.
*
* @param  divisor Value by which this {@code Complex} is to be divided.
* @return {@code this / divisor}.
* @see #divide(Complex)
*/
Complex Complex::operator/(double divisor)
//Complex Complex::divide(double divisor)
{
	if (IsNaN || isnan(divisor))
	{
		return NaN;
	}
	if (divisor == 0)
	{
		return NaN;
	}
	if (isinf(divisor))
	{
		return !IsInfinite ? ZERO : NaN;
	}
	return createComplex(real / divisor,
		imaginary / divisor);
}

/** {@inheritDoc} */
Complex Complex::reciprocal()
{
	if (IsNaN)
	{
		return NaN;
	}

	if (real == 0.0 && imaginary == 0.0)
	{
		return INF;
	}

	if (IsInfinite)
	{
		return ZERO;
	}

	if (fabs(real) < fabs(imaginary))
	{
		double q = real / imaginary;
		double scale = 1.0 / (real * q + imaginary);
		return createComplex(scale * q, -scale);
	}
	else
	{
		double q = imaginary / real;
		double scale = 1.0 / (imaginary * q + real);
		return createComplex(scale, -scale * q);
	}
}

/**
* Test for equality with another object.
* If both the real and imaginary parts of two Complex numbers
* are exactly the same, and neither is {@code Double.NaN}, the two
* Complex objects are considered to be equal.
* The behavior is the same as for JDK's {@link Double#equals(Object)
* Double}:
* <ul>
*  <li>All {@code NaN} values are considered to be equal,
*   i.e, if either (or both) real and imaginary parts of the Complex
*   number are equal to {@code Double.NaN}, the Complex number is equal
*   to {@code NaN}.
*  </li>
*  <li>
*   Instances constructed with different representations of zero (i.e.
*   either "0" or "-0") are <em>not</em> considered to be equal.
*  </li>
* </ul>
*
* @param other Object to test for equality with this instance.
* @return {@code true} if the objects are equal, {@code false} if object
* is {@code null}, not an instance of {@code Complex}, or not equal to
* this instance.
*/
bool Complex::Equals(const Complex &c) const
{
	if (c.IsNaN)
	{
		return IsNaN;
	}
	else
	{
		return real == c.real &&
			imaginary == c.imaginary;
	}
}

/**
* Returns {@code true} if, both for the real part and for the imaginary
* part, there is no double value strictly between the arguments or the
* difference between them is within the range of allowed error
* (inclusive).  Returns {@code false} if either of the arguments is NaN.
*
* @param x First value (cannot be {@code null}).
* @param y Second value (cannot be {@code null}).
* @param eps Amount of allowed absolute error.
* @return {@code true} if the values are two adjacent floating point
* numbers or they are within range of each other.
*
* @see Precision#equals(double,double,double)
*/
bool Complex::equals(const Complex &x, const Complex &y, double eps)
{
	return fabs(x.real - y.real) <= eps &&
		fabs(x.imaginary - y.imaginary) <= eps;
}

/**
* Returns a {@code Complex} whose value is {@code this * factor}.
* Implements preliminary checks for {@code NaN} and infinity followed by
* the definitional formula:
* <p>
*   {@code (a + bi)(c + di) = (ac - bd) + (ad + bc)i}
* </p>
* Returns {@link #NaN} if either {@code this} or {@code factor} has one or
* more {@code NaN} parts.
* <p>
* Returns {@link #INF} if neither {@code this} nor {@code factor} has one
* or more {@code NaN} parts and if either {@code this} or {@code factor}
* has one or more infinite parts (same result is returned regardless of
* the sign of the components).
* </p><p>
* Returns finite values in components of the result per the definitional
* formula in all remaining cases.</p>
*
* @param  factor value to be multiplied by this {@code Complex}.
* @return {@code this * factor}.
*/
Complex Complex::operator*(const Complex &factor)
{
	if (IsNaN || factor.IsNaN)
	{
		return NaN;
	}
	if (isinf(real) ||
		isinf(imaginary) ||
		isinf(factor.real) ||
		isinf(factor.imaginary))
	{
		// we don't use IsInfinite to avoid testing for NaN again
		return INF;
	}
	return createComplex(real * factor.real - imaginary * factor.imaginary,
		real * factor.imaginary + imaginary * factor.real);
}

/**
* Returns a {@code Complex} whose value is {@code this * factor}, with {@code factor}
* interpreted as a integer number.
*
* @param  factor value to be multiplied by this {@code Complex}.
* @return {@code this * factor}.
* @see #multiply(Complex)
*/
Complex Complex::operator*(int factor)
//Complex Complex::multiply(int factor)
{
	if (IsNaN)
	{
		return NaN;
	}
	if (isinf(real) ||
		isinf(imaginary))
	{
		return INF;
	}
	return createComplex(real * factor, imaginary * factor);
}

/**
* Returns a {@code Complex} whose value is {@code this * factor}, with {@code factor}
* interpreted as a real number.
*
* @param  factor value to be multiplied by this {@code Complex}.
* @return {@code this * factor}.
* @see #multiply(Complex)
*/
Complex Complex::operator*(double factor)
//Complex Complex::multiply(double factor)
{
	if (IsNaN || isnan(factor))
	{
		return NaN;
	}
	if (isinf(real) ||
		isinf(imaginary) ||
		isinf(factor))
	{
		// we don't use IsInfinite to avoid testing for NaN again
		return INF;
	}
	return createComplex(real * factor, imaginary * factor);
}

/**
* Returns a {@code Complex} whose value is {@code (-this)}.
* Returns {@code NaN} if either real or imaginary
* part of this Complex number is {@code Double.NaN}.
*
* @return {@code -this}.
*/
Complex Complex::operator-()
//Complex Complex::negate()
{
	if (IsNaN)
	{
		return NaN;
	}
	return createComplex(-real, -imaginary);
}

/**
* Returns a {@code Complex} whose value is
* {@code (this - subtrahend)}.
* Uses the definitional formula
* <p>
*  {@code (a + bi) - (c + di) = (a-c) + (b-d)i}
* </p>
* If either {@code this} or {@code subtrahend} has a {@code NaN]} value in either part,
* {@link #NaN} is returned; otherwise infinite and {@code NaN} values are
* returned in the parts of the result according to the rules for
* {@link java.lang.Double} arithmetic.
*
* @param  subtrahend value to be subtracted from this {@code Complex}.
* @return {@code this - subtrahend}.
*/
Complex Complex::operator-(const Complex &subtrahend)
{
	if (IsNaN || subtrahend.IsNaN)
	{
		return NaN;
	}
	return createComplex(real - subtrahend.getReal(),
		imaginary - subtrahend.getImaginary());
}

/**
* Returns a {@code Complex} whose value is
* {@code (this - subtrahend)}.
*
* @param  subtrahend value to be subtracted from this {@code Complex}.
* @return {@code this - subtrahend}.
* @see #subtract(Complex)
*/
Complex Complex::operator-(double subtrahend)
//Complex Complex::subtract(double subtrahend)
{
	if (IsNaN || isnan(subtrahend))
	{
		return NaN;
	}
	return createComplex(real - subtrahend, imaginary);
}

/**
* Compute the
* <a href="http://mathworld.wolfram.com/InverseCosine.html" TARGET="_top">
* inverse cosine</a> of this Complex number.
* Implements the formula:
* <p>
*  {@code acos(z) = -i (log(z + i (sqrt(1 - z<sup>2</sup>))))}
* </p>
* Returns {@link Complex#NaN} if either real or imaginary part of the
* input argument is {@code NaN} or infinite.
*
* @return the inverse cosine of this Complex number.
*/
Complex Complex::acos()
{
	if (IsNaN)
	{
		return NaN;
	}
	return - (((*this) + (*this).sqrt1z() * I).log())*I;
}

/**
* Compute the
* <a href="http://mathworld.wolfram.com/InverseSine.html" TARGET="_top">
* inverse sine</a> of this Complex number.
* Implements the formula:
* <p>
*  {@code asin(z) = -i (log(sqrt(1 - z<sup>2</sup>) + iz))}
* </p><p>
* Returns {@link Complex#NaN} if either real or imaginary part of the
* input argument is {@code NaN} or infinite.</p>
*
* @return the inverse sine of this Complex number.
*/
Complex Complex::asin()
{
	if (IsNaN)
	{
		return NaN;
	}
	return (this->sqrt1z() + ((*this)*I)).log() * (-I);
}

/**
* Compute the
* <a href="http://mathworld.wolfram.com/InverseTangent.html" TARGET="_top">
* inverse tangent</a> of this Complex number.
* Implements the formula:
* <p>
* {@code atan(z) = (i/2) log((i + z)/(i - z))}
* </p><p>
* Returns {@link Complex#NaN} if either real or imaginary part of the
* input argument is {@code NaN} or infinite.</p>
*
* @return the inverse tangent of this Complex number
*/
Complex Complex::atan()
{
	if (IsNaN)
	{
		return NaN;
	}
	return ((I + (*this)) / (I - (*this))).log() * I *0.5;
}

/**
* Compute the
* <a href="http://mathworld.wolfram.com/Cosine.html" TARGET="_top">
* cosine</a> of this Complex number.
* Implements the formula:
* <p>
*  {@code cos(a + bi) = cos(a)cosh(b) - sin(a)sinh(b)i}
* </p><p>
* where the (real) functions on the right-hand side are
* {@link FastMath#sin}, {@link FastMath#cos},
* {@link FastMath#cosh} and {@link FastMath#sinh}.
* </p><p>
* Returns {@link Complex#NaN} if either real or imaginary part of the
* input argument is {@code NaN}.
* </p><p>
* Infinite values in real or imaginary parts of the input may result in
* infinite or NaN values returned in parts of the result.</p>
* <pre>
*  Examples:
*  <code>
*   cos(1 &plusmn; INFINITY i) = 1 \u2213 INFINITY i
*   cos(&plusmn;INFINITY + i) = NaN + NaN i
*   cos(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
*  </code>
* </pre>
*
* @return the cosine of this Complex number.
*/
Complex Complex::cos()
{
	if (IsNaN)
	{
		return NaN;
	}
	return createComplex(::cos(real) * ::cosh(imaginary),
		-::sin(real) * ::sinh(imaginary));
}

/**
* Compute the
* <a href="http://mathworld.wolfram.com/HyperbolicCosine.html" TARGET="_top">
* hyperbolic cosine</a> of this Complex number.
* Implements the formula:
* <pre>
*  <code>
*   cosh(a + bi) = cosh(a)cos(b) + sinh(a)sin(b)i
*  </code>
* </pre>
* where the (real) functions on the right-hand side are
* {@link FastMath#sin}, {@link FastMath#cos},
* {@link FastMath#cosh} and {@link FastMath#sinh}.
* <p>
* Returns {@link Complex#NaN} if either real or imaginary part of the
* input argument is {@code NaN}.
* </p>
* Infinite values in real or imaginary parts of the input may result in
* infinite or NaN values returned in parts of the result.
* <pre>
*  Examples:
*  <code>
*   cosh(1 &plusmn; INFINITY i) = NaN + NaN i
*   cosh(&plusmn;INFINITY + i) = INFINITY &plusmn; INFINITY i
*   cosh(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
*  </code>
* </pre>
*
* @return the hyperbolic cosine of this Complex number.
*/
Complex Complex::cosh()
{
	if (IsNaN)
	{
		return NaN;
	}
	return createComplex(::cosh(real) * ::cos(imaginary),
		::sinh(real) * ::sin(imaginary));
}

/**
* Compute the
* <a href="http://mathworld.wolfram.com/ExponentialFunction.html" TARGET="_top">
* exponential function</a> of this Complex number.
* Implements the formula:
* <pre>
*  <code>
*   exp(a + bi) = exp(a)cos(b) + exp(a)sin(b)i
*  </code>
* </pre>
* where the (real) functions on the right-hand side are
* {@link FastMath#exp}, {@link FastMath#cos}, and
* {@link FastMath#sin}.
* <p>
* Returns {@link Complex#NaN} if either real or imaginary part of the
* input argument is {@code NaN}.
* </p>
* Infinite values in real or imaginary parts of the input may result in
* infinite or NaN values returned in parts of the result.
* <pre>
*  Examples:
*  <code>
*   exp(1 &plusmn; INFINITY i) = NaN + NaN i
*   exp(INFINITY + i) = INFINITY + INFINITY i
*   exp(-INFINITY + i) = 0 + 0i
*   exp(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
*  </code>
* </pre>
*
* @return <code><i>e</i><sup>this</sup></code>.
*/
Complex Complex::exp()
{
	if (IsNaN)
	{
		return NaN;
	}
	double expReal = ::exp(real);
	return createComplex(expReal * ::cos(imaginary),
		expReal * ::sin(imaginary));
}

/**
* Compute the
* <a href="http://mathworld.wolfram.com/NaturalLogarithm.html" TARGET="_top">
* natural logarithm</a> of this Complex number.
* Implements the formula:
* <pre>
*  <code>
*   log(a + bi) = ln(|a + bi|) + arg(a + bi)i
*  </code>
* </pre>
* where ln on the right hand side is {@link FastMath#log},
* {@code |a + bi|} is the modulus, {@link Complex#abs},  and
* {@code arg(a + bi) = }{@link FastMath#atan2}(b, a).
* <p>
* Returns {@link Complex#NaN} if either real or imaginary part of the
* input argument is {@code NaN}.
* </p>
* Infinite (or critical) values in real or imaginary parts of the input may
* result in infinite or NaN values returned in parts of the result.
* <pre>
*  Examples:
*  <code>
*   log(1 &plusmn; INFINITY i) = INFINITY &plusmn; (&pi;/2)i
*   log(INFINITY + i) = INFINITY + 0i
*   log(-INFINITY + i) = INFINITY + &pi;i
*   log(INFINITY &plusmn; INFINITY i) = INFINITY &plusmn; (&pi;/4)i
*   log(-INFINITY &plusmn; INFINITY i) = INFINITY &plusmn; (3&pi;/4)i
*   log(0 + 0i) = -INFINITY + 0i
*  </code>
* </pre>
*
* @return the value <code>ln &nbsp; this</code>, the natural logarithm
* of {@code this}.
*/
Complex Complex::log()
{
	if (IsNaN)
	{
		return NaN;
	}
	return createComplex(::log(abs()),
		::atan2(imaginary, real));
}

/**
* Returns of value of this Complex number raised to the power of {@code x}.
* Implements the formula:
* <pre>
*  <code>
*   y<sup>x</sup> = exp(x&middot;log(y))
*  </code>
* </pre>
* where {@code exp} and {@code log} are {@link #exp} and
* {@link #log}, respectively.
* <p>
* Returns {@link Complex#NaN} if either real or imaginary part of the
* input argument is {@code NaN} or infinite, or if {@code y}
* equals {@link Complex#ZERO}.</p>
*
* @param  x exponent to which this {@code Complex} is to be raised.
* @return <code> this<sup>x</sup></code>.     *
*/
Complex Complex::pow(const Complex &x)
{
	return (this->log() * x).exp();
}

/**
* Returns of value of this Complex number raised to the power of {@code x}.
*
* @param  x exponent to which this {@code Complex} is to be raised.
* @return <code>this<sup>x</sup></code>.
* @see #pow(Complex)
*/
Complex Complex::pow(double x)
{
	return (this->log()* x).exp();
}

/**
* Compute the
* <a href="http://mathworld.wolfram.com/Sine.html" TARGET="_top">
* sine</a>
* of this Complex number.
* Implements the formula:
* <pre>
*  <code>
*   sin(a + bi) = sin(a)cosh(b) - cos(a)sinh(b)i
*  </code>
* </pre>
* where the (real) functions on the right-hand side are
* {@link FastMath#sin}, {@link FastMath#cos},
* {@link FastMath#cosh} and {@link FastMath#sinh}.
* <p>
* Returns {@link Complex#NaN} if either real or imaginary part of the
* input argument is {@code NaN}.
* </p><p>
* Infinite values in real or imaginary parts of the input may result in
* infinite or {@code NaN} values returned in parts of the result.
* <pre>
*  Examples:
*  <code>
*   sin(1 &plusmn; INFINITY i) = 1 &plusmn; INFINITY i
*   sin(&plusmn;INFINITY + i) = NaN + NaN i
*   sin(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
*  </code>
* </pre>
*
* @return the sine of this Complex number.
*/
Complex Complex::sin()
{
	if (IsNaN)
	{
		return NaN;
	}
	return createComplex(::sin(real) * ::cosh(imaginary),
		::cos(real) * ::sinh(imaginary));
}

/**
* Compute the
* <a href="http://mathworld.wolfram.com/HyperbolicSine.html" TARGET="_top">
* hyperbolic sine</a> of this Complex number.
* Implements the formula:
* <pre>
*  <code>
*   sinh(a + bi) = sinh(a)cos(b)) + cosh(a)sin(b)i
*  </code>
* </pre>
* where the (real) functions on the right-hand side are
* {@link FastMath#sin}, {@link FastMath#cos},
* {@link FastMath#cosh} and {@link FastMath#sinh}.
* <p>
* Returns {@link Complex#NaN} if either real or imaginary part of the
* input argument is {@code NaN}.
* </p><p>
* Infinite values in real or imaginary parts of the input may result in
* infinite or NaN values returned in parts of the result.
* <pre>
*  Examples:
*  <code>
*   sinh(1 &plusmn; INFINITY i) = NaN + NaN i
*   sinh(&plusmn;INFINITY + i) = &plusmn; INFINITY + INFINITY i
*   sinh(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
*  </code>
* </pre>
*
* @return the hyperbolic sine of {@code this}.
*/
Complex Complex::sinh()
{
	if (IsNaN)
	{
		return NaN;
	}
	return createComplex(::sinh(real) * ::cos(imaginary),
		::cosh(real) * ::sin(imaginary));
}

/**
* Compute the
* <a href="http://mathworld.wolfram.com/SquareRoot.html" TARGET="_top">
* square root</a> of this Complex number.
* Implements the following algorithm to compute {@code sqrt(a + bi)}:
* <ol><li>Let {@code t = sqrt((|a| + |a + bi|) / 2)}</li>
* <li><pre>if {@code  a &#8805; 0} return {@code t + (b/2t)i}
*  else return {@code |b|/2t + sign(b)t i }</pre></li>
* </ol>
* where <ul>
* <li>{@code |a| = }{@link FastMath#abs}(a)</li>
* <li>{@code |a + bi| = }{@link Complex#abs}(a + bi)</li>
* <li>{@code sign(b) =  }{@link FastMath#copySign(double,double) copySign(1d, b)}
* </ul>
* <p>
* Returns {@link Complex#NaN} if either real or imaginary part of the
* input argument is {@code NaN}.
* </p>
* Infinite values in real or imaginary parts of the input may result in
* infinite or NaN values returned in parts of the result.
* <pre>
*  Examples:
*  <code>
*   sqrt(1 &plusmn; INFINITY i) = INFINITY + NaN i
*   sqrt(INFINITY + i) = INFINITY + 0i
*   sqrt(-INFINITY + i) = 0 + INFINITY i
*   sqrt(INFINITY &plusmn; INFINITY i) = INFINITY + NaN i
*   sqrt(-INFINITY &plusmn; INFINITY i) = NaN &plusmn; INFINITY i
*  </code>
* </pre>
*
* @return the square root of {@code this}.
*/
Complex Complex::sqrt()
{
	if (IsNaN)
	{
		return NaN;
	}

	if (real == 0.0 && imaginary == 0.0)
	{
		return createComplex(0.0, 0.0);
	}

	double t = ::sqrt((fabs(real) + abs()) / 2.0);
	if (real >= 0.0)
	{
		return createComplex(t, imaginary / (2.0 * t));
	}
	else
	{
		double signedT = imaginary < 0 ? -t : t;
		return createComplex(fabs(imaginary) / (2.0 * t),
			signedT);
	}
}

/**
* Compute the
* <a href="http://mathworld.wolfram.com/SquareRoot.html" TARGET="_top">
* square root</a> of <code>1 - this<sup>2</sup></code> for this Complex
* number.
* Computes the result directly as
* {@code sqrt(ONE.subtract(z.multiply(z)))}.
* <p>
* Returns {@link Complex#NaN} if either real or imaginary part of the
* input argument is {@code NaN}.
* </p>
* Infinite values in real or imaginary parts of the input may result in
* infinite or NaN values returned in parts of the result.
*
* @return the square root of <code>1 - this<sup>2</sup></code>.
*/
Complex Complex::sqrt1z()
{
	return (ONE - ((*this) * (*this))).sqrt();
}

/**
* Compute the
* <a href="http://mathworld.wolfram.com/Tangent.html" TARGET="_top">
* tangent</a> of this Complex number.
* Implements the formula:
* <pre>
*  <code>
*   tan(a + bi) = sin(2a)/(cos(2a)+cosh(2b)) + [sinh(2b)/(cos(2a)+cosh(2b))]i
*  </code>
* </pre>
* where the (real) functions on the right-hand side are
* {@link FastMath#sin}, {@link FastMath#cos}, {@link FastMath#cosh} and
* {@link FastMath#sinh}.
* <p>
* Returns {@link Complex#NaN} if either real or imaginary part of the
* input argument is {@code NaN}.
* </p>
* Infinite (or critical) values in real or imaginary parts of the input may
* result in infinite or NaN values returned in parts of the result.
* <pre>
*  Examples:
*  <code>
*   tan(a &plusmn; INFINITY i) = 0 &plusmn; i
*   tan(&plusmn;INFINITY + bi) = NaN + NaN i
*   tan(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
*   tan(&plusmn;&pi;/2 + 0 i) = &plusmn;INFINITY + NaN i
*  </code>
* </pre>
*
* @return the tangent of {@code this}.
*/
Complex Complex::tan()
{
	if (IsNaN || isinf(real))
	{
		return NaN;
	}
	if (imaginary > 20.0)
	{
		return createComplex(0.0, 1.0);
	}
	if (imaginary < -20.0)
	{
		return createComplex(0.0, -1.0);
	}

	double real2 = 2.0 * real;
	double imaginary2 = 2.0 * imaginary;
	double d = ::cos(real2) + ::cosh(imaginary2);

	return createComplex(::sin(real2) / d,
		::sinh(imaginary2) / d);
}

/**
* Compute the
* <a href="http://mathworld.wolfram.com/HyperbolicTangent.html" TARGET="_top">
* hyperbolic tangent</a> of this Complex number.
* Implements the formula:
* <pre>
*  <code>
*   tan(a + bi) = sinh(2a)/(cosh(2a)+cos(2b)) + [sin(2b)/(cosh(2a)+cos(2b))]i
*  </code>
* </pre>
* where the (real) functions on the right-hand side are
* {@link FastMath#sin}, {@link FastMath#cos}, {@link FastMath#cosh} and
* {@link FastMath#sinh}.
* <p>
* Returns {@link Complex#NaN} if either real or imaginary part of the
* input argument is {@code NaN}.
* </p>
* Infinite values in real or imaginary parts of the input may result in
* infinite or NaN values returned in parts of the result.
* <pre>
*  Examples:
*  <code>
*   tanh(a &plusmn; INFINITY i) = NaN + NaN i
*   tanh(&plusmn;INFINITY + bi) = &plusmn;1 + 0 i
*   tanh(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
*   tanh(0 + (&pi;/2)i) = NaN + INFINITY i
*  </code>
* </pre>
*
* @return the hyperbolic tangent of {@code this}.
*/
Complex Complex::tanh()
{
	if (IsNaN || isinf(imaginary))
	{
		return NaN;
	}
	if (real > 20.0)
	{
		return createComplex(1.0, 0.0);
	}
	if (real < -20.0)
	{
		return createComplex(-1.0, 0.0);
	}
	double real2 = 2.0 * real;
	double imaginary2 = 2.0 * imaginary;
	double d = ::cosh(real2) + ::cos(imaginary2);

	return createComplex(::sinh(real2) / d,
		::sin(imaginary2) / d);
}


/**
* Compute the argument of this Complex number.
* The argument is the angle phi between the positive real axis and
* the point representing this number in the Complex plane.
* The value returned is between -PI (not inclusive)
* and PI (inclusive), with negative values returned for numbers with
* negative imaginary parts.
* <p>
* If either real or imaginary part (or both) is NaN, NaN is returned.
* Infinite parts are handled as {@code Math.atan2} handles them,
* essentially treating finite parts as zero in the presence of an
* infinite coordinate and returning a multiple of pi/4 depending on
* the signs of the infinite parts.
* See the javadoc for {@code Math.atan2} for full details.
*
* @return the argument of {@code this}.
*/
double Complex::getArgument()
{
	return ::atan2(getImaginary(), getReal());
}

/**
* Computes first the n-th root of this Complex number.
* The nth roots are defined by the formula:
* <pre>
*  <code>
*   z<sub>k</sub> = abs<sup>1/n</sup> (cos(phi + 2&pi;k/n) + i (sin(phi + 2&pi;k/n))
*  </code>
* </pre>
* for <i>{@code k=0, 1, ..., n-1}</i>, where {@code abs} and {@code phi}
* are respectively the {@link #abs() modulus} and
* {@link #getArgument() argument} of this Complex number.
* <p>
* If one or both parts of this Complex number is NaN, a list with just
* one element, {@link #NaN} is returned.
* if neither part is NaN, but at least one part is infinite, the result
* is a one-element list containing {@link #INF}.
*
* @param n Degree of root.
* @return first {@code n}-th roots of {@code this}.     *
*/
Complex Complex::nthRoot(int n)
{
	if (n <= 0)
	{
		throw "Can not compute nth root for negative n";
	}
	if (IsNaN) return NaN;
	if (IsInfinite) return INF;

	// nth root of abs -- faster / more accurate to use a solver here?
	double nthRootOfAbs = ::pow(abs(), 1.0 / n);

	// Compute nth roots of Complex number with k = 0, 1, ... n-1
	double nthPhi = getArgument() / n;
	double innerPart = nthPhi;
	double realPart = nthRootOfAbs * ::cos(innerPart);
	double imaginaryPart = nthRootOfAbs * ::sin(innerPart);
	return createComplex(realPart, imaginaryPart);
}

/**
* Computes the n-th roots of this Complex number.
* The nth roots are defined by the formula:
* <pre>
*  <code>
*   z<sub>k</sub> = abs<sup>1/n</sup> (cos(phi + 2&pi;k/n) + i (sin(phi + 2&pi;k/n))
*  </code>
* </pre>
* for <i>{@code k=0, 1, ..., n-1}</i>, where {@code abs} and {@code phi}
* are respectively the {@link #abs() modulus} and
* {@link #getArgument() argument} of this Complex number.
* <p>
* If one or both parts of this Complex number is NaN, a list with just
* one element, {@link #NaN} is returned.
* if neither part is NaN, but at least one part is infinite, the result
* is a one-element list containing {@link #INF}.
*
* @param n Degree of root.
* @return array all {@code n}-th roots of {@code this}.     *
*/
vector<Complex> Complex::nthRoots(int n)
{
	if (n <= 0)
	{
		throw "Can not compute nth root for negative n";
	}
	vector<Complex> result;

	if (IsNaN)
	{
		for (int k = 0; k < n; k++) result[k] = NaN;
		return result;
	}
	if (IsInfinite)
	{
		for (int k = 0; k < n; k++) result[k] = INF;
		return result;
	}

	// nth root of abs -- faster / more accurate to use a solver here?
	double nthRootOfAbs = ::pow(abs(), 1.0 / n);

	// Compute nth roots of Complex number with k = 0, 1, ... n-1
	double nthPhi = getArgument() / n;
	double slice = 2 * M_PI / n;
	double innerPart = nthPhi;
	for (int k = 0; k < n; k++)
	{
		// inner part
		double realPart = nthRootOfAbs * ::cos(innerPart);
		double imaginaryPart = nthRootOfAbs * ::sin(innerPart);
		result.push_back(createComplex(realPart, imaginaryPart));
		innerPart += slice;
	}
	return result;
}

/**
* Create a Complex number given the real and imaginary parts.
*
* @param realPart Real part.
* @param imaginaryPart Imaginary part.
* @return a new Complex number instance.
* @see #valueOf(double, double)
*/
Complex Complex::createComplex(const double realPart,const double imaginaryPart)
{
	Complex c(realPart, imaginaryPart);
	return c;
}

string Complex::toString()
{
	return "(" + to_string(real) + ", " + to_string(imaginary) + ")";
}

using System;

namespace Util
{
    class Complex
    {
        /** The square root of -1. A number representing "0.0 + 1.0i" */
        public static Complex I = new Complex(0.0, 1.0);
        // CHECKSTYLE: stop ConstantName
        /** A Complex number representing "NaN + NaNi" */
        public static Complex NaN = new Complex(Double.NaN, Double.NaN);
        // CHECKSTYLE: resume ConstantName
        /** A Complex number representing "+INF + INFi" */
        public static Complex INF = new Complex(Double.PositiveInfinity, Double.PositiveInfinity);
        /** A Complex number representing "1.0 + 0.0i" */
        public static Complex ONE = new Complex(1.0, 0.0);
        /** A Complex number representing "0.0 + 0.0i" */
        public static Complex ZERO = new Complex(0.0, 0.0);

        /** The imaginary part. */
        private double imaginary;
        /** The real part. */
        private double real;
        /** Record whether this Complex number is equal to NaN. */
        private bool IsNaN;
        /** Record whether this Complex number is infinite. */
        private bool IsInfinite;

        /**
         * Create a Complex number given only the real part.
         *
         * @param real Real part.
         */
        public Complex(double real) : this(real, 0.0)
        { }

        /**
         * Create a Complex number given the real and imaginary parts.
         *
         * @param real Real part.
         * @param imaginary Imaginary part.
         */
        public Complex(double real, double imaginary)
        {
            this.real = real;
            this.imaginary = imaginary;

            IsNaN = Double.IsNaN(real) || Double.IsNaN(imaginary);
            IsInfinite = !IsNaN &&
                (Double.IsInfinity(real) || Double.IsInfinity(imaginary));
        }

        /**
         * Return the absolute value of this Complex number.
         * Returns {@code NaN} if either real or imaginary part is {@code NaN}
         * and {@code Double.POSITIVE_INFINITY} if neither part is {@code NaN},
         * but at least one part is infinite.
         *
         * @return the absolute value.
         */
        public double abs()
        {
            if (IsNaN)
            {
                return Double.NaN;
            }
            if (IsInfinite)
            {
                return Double.PositiveInfinity;
            }
            if (Math.Abs(real) < Math.Abs(imaginary))
            {
                if (imaginary == 0.0)
                {
                    return Math.Abs(real);
                }
                double q = real / imaginary;
                return Math.Abs(imaginary) * Math.Sqrt(1 + q * q);
            }
            else
            {
                if (real == 0.0)
                {
                    return Math.Abs(imaginary);
                }
                double q = imaginary / real;
                return Math.Abs(real) * Math.Sqrt(1 + q * q);
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
        public static Complex operator+(Complex first, Complex addend)
        {
            return first.add(addend);
        }

        public Complex add(Complex addend)
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
        public static Complex operator +(Complex first, double addend)
        {
            return first.add(addend);
        }

        public Complex add(double addend)
        {
            if (IsNaN || Double.IsNaN(addend))
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
        public Complex conjugate()
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
        public static Complex operator /(Complex first, Complex divisor)
        {
            return first.divide(divisor);
        }

        public Complex divide(Complex divisor)
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

            if (Math.Abs(c) < Math.Abs(d))
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
        public static Complex operator /(Complex first, double divisor)
        {
            return first.divide(divisor);
        }

        public Complex divide(double divisor)
        {
            if (IsNaN || Double.IsNaN(divisor))
            {
                return NaN;
            }
            if (divisor == 0d)
            {
                return NaN;
            }
            if (Double.IsInfinity(divisor))
            {
                return !IsInfinite ? ZERO : NaN;
            }
            return createComplex(real / divisor,
                                 imaginary / divisor);
        }

        /** {@inheritDoc} */
        public Complex reciprocal()
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

            if (Math.Abs(real) < Math.Abs(imaginary))
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
    public override bool Equals(Object other)
        {
            if (this == other)
            {
                return true;
            }
            if (other is Complex){
                Complex c = (Complex)other;
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
            return false;
        }

        /**
         * Returns {@code true} iff the values are equal as defined by
         * {@link #equals(Complex,Complex,int) equals(x, y, 1)}.
         *
         * @param x First value (cannot be {@code null}).
         * @param y Second value (cannot be {@code null}).
         * @return {@code true} if the values are equal.
         */
        public static bool equals(Complex x, Complex y)
        {
            return equals(x, y, 1);
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
        public static bool equals(Complex x, Complex y, double eps)
        {
            return Math.Abs(x.real - y.real) <= eps &&
                    Math.Abs(x.imaginary - y.imaginary) <= eps;
        }

        /**
         * Access the imaginary part.
         *
         * @return the imaginary part.
         */
        public double getImaginary()
        {
            return imaginary;
        }

        /**
         * Access the real part.
         *
         * @return the real part.
         */
        public double getReal()
        {
            return real;
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
        public static Complex operator *(Complex first, Complex factor)
        {
            return first.multiply(factor);
        }

        public Complex multiply(Complex factor)
        {
            if (IsNaN || factor.IsNaN)
            {
                return NaN;
            }
            if (Double.IsInfinity(real) ||
                Double.IsInfinity(imaginary) ||
                Double.IsInfinity(factor.real) ||
                Double.IsInfinity(factor.imaginary))
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
        public static Complex operator*(Complex first, int factor)
        {
            return first.multiply(factor);
        }

        public Complex multiply(int factor)
        {
            if (IsNaN)
            {
                return NaN;
            }
            if (Double.IsInfinity(real) ||
                Double.IsInfinity(imaginary))
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
        public static Complex operator *(Complex first, double factor)
        {
            return first.multiply(factor);
        }

        public Complex multiply(double factor)
        {
            if (IsNaN || Double.IsNaN(factor))
            {
                return NaN;
            }
            if (Double.IsInfinity(real) ||
                Double.IsInfinity(imaginary) ||
                Double.IsInfinity(factor))
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
        public static Complex operator -(Complex first)
        {
            return first.negate();
        }

        public Complex negate()
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
        public static Complex operator -(Complex first, Complex subtrahend)
        {
            return first.subtract(subtrahend);
        }

        public Complex subtract(Complex subtrahend)
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
        public static Complex operator -(Complex first, double subtrahend)
        {
            return first.subtract(subtrahend);
        }

        public Complex subtract(double subtrahend)
        {
            if (IsNaN || Double.IsNaN(subtrahend))
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
        public Complex acos()
        {
            if (IsNaN)
            {
                return NaN;
            }

            return this.add(this.sqrt1z().multiply(I)).log().multiply(I.negate());
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
        public Complex asin()
        {
            if (IsNaN)
            {
                return NaN;
            }
            return this.sqrt1z().add(this.multiply(I)).log().multiply(I.negate());
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
        public Complex atan()
        {
            if (IsNaN)
            {
                return NaN;
            }

            return this.add(I).divide(I.subtract(this)).log()
                    .multiply(I.divide(createComplex(2.0, 0.0)));
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
        public Complex cos()
        {
            if (IsNaN)
            {
                return NaN;
            }
            return createComplex(Math.Cos(real) * Math.Cosh(imaginary),
                                 -Math.Sin(real) * Math.Sinh(imaginary));
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
        public Complex cosh()
        {
            if (IsNaN)
            {
                return NaN;
            }
            return createComplex(Math.Cosh(real) * Math.Cos(imaginary),
                                 Math.Sinh(real) * Math.Sin(imaginary));
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
        public Complex exp()
        {
            if (IsNaN)
            {
                return NaN;
            }

            double expReal = Math.Exp(real);
            return createComplex(expReal * Math.Cos(imaginary),
                                 expReal * Math.Sin(imaginary));
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
        public Complex log()
        {
            if (IsNaN)
            {
                return NaN;
            }
            return createComplex(Math.Log(abs()),
                                 Math.Atan2(imaginary, real));
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
        public Complex pow(Complex x)
        {
            return this.log().multiply(x).exp();
        }

        /**
         * Returns of value of this Complex number raised to the power of {@code x}.
         *
         * @param  x exponent to which this {@code Complex} is to be raised.
         * @return <code>this<sup>x</sup></code>.
         * @see #pow(Complex)
         */
        public Complex pow(double x)
        {
            return this.log().multiply(x).exp();
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
        public Complex sin()
        {
            if (IsNaN)
            {
                return NaN;
            }
            return createComplex(Math.Sin(real) * Math.Cosh(imaginary),
                                 Math.Cos(real) * Math.Sinh(imaginary));
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
        public Complex sinh()
        {
            if (IsNaN)
            {
                return NaN;
            }
            return createComplex(Math.Sinh(real) * Math.Cos(imaginary),
                Math.Cosh(real) * Math.Sin(imaginary));
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
        public Complex sqrt()
        {
            if (IsNaN)
            {
                return NaN;
            }

            if (real == 0.0 && imaginary == 0.0)
            {
                return createComplex(0.0, 0.0);
            }

            double t = Math.Sqrt((Math.Abs(real) + abs()) / 2.0);
            if (real >= 0.0)
            {
                return createComplex(t, imaginary / (2.0 * t));
            }
            else
            {
                double signedT = imaginary < 0 ? -t : t;
                return createComplex(Math.Abs(imaginary) / (2.0 * t),
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
        public Complex sqrt1z()
        {
            return createComplex(1.0, 0.0).subtract(this.multiply(this)).sqrt();
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
        public Complex tan()
        {
            if (IsNaN || Double.IsInfinity(real))
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
            double d = Math.Cos(real2) + Math.Cosh(imaginary2);

            return createComplex(Math.Sin(real2) / d,
                                 Math.Sinh(imaginary2) / d);
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
        public Complex tanh()
        {
            if (IsNaN || Double.IsInfinity(imaginary))
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
            double d = Math.Cosh(real2) + Math.Cos(imaginary2);

            return createComplex(Math.Sinh(real2) / d,
                                 Math.Sin(imaginary2) / d);
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
        public double getArgument()
        {
            return Math.Atan2(getImaginary(), getReal());
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
        public Complex nthRoot(int n)
        {
            if (n <= 0)
            {
                throw new Exception("Can not compute nth root for negative n");
            }
            if (IsNaN) return NaN;
            if (IsInfinite) return INF;

            // nth root of abs -- faster / more accurate to use a solver here?
            double nthRootOfAbs = Math.Pow(abs(), 1.0 / n);

            // Compute nth roots of Complex number with k = 0, 1, ... n-1
            double nthPhi = getArgument() / n;
            double innerPart = nthPhi;
            double realPart = nthRootOfAbs * Math.Cos(innerPart);
            double imaginaryPart = nthRootOfAbs * Math.Sin(innerPart);
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
        public Complex[] nthRoots(int n)
        {

            if (n <= 0)
            {
                throw new Exception("Can not compute nth root for negative n");
            }

            Complex[]
            result = new Complex[n];

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
            double nthRootOfAbs = Math.Pow(abs(), 1.0 / n);

            // Compute nth roots of Complex number with k = 0, 1, ... n-1
            double nthPhi = getArgument() / n;
            double slice = 2 * Math.PI / n;
            double innerPart = nthPhi;
            for (int k = 0; k < n; k++)
            {
                // inner part
                double realPart = nthRootOfAbs * Math.Cos(innerPart);
                double imaginaryPart = nthRootOfAbs * Math.Sin(innerPart);
                result[k] = createComplex(realPart, imaginaryPart);
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
        public static Complex createComplex(double realPart,
                                        double imaginaryPart)
        {
            return new Complex(realPart, imaginaryPart);
        }

        public override String ToString()
        {
            return "(" + real + ", " + imaginary + ")";
        }
    }
}

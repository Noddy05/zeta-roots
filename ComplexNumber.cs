using System;

namespace ComplexNumbers
{
    struct ComplexNumber
    {
        #region Instantiation
        public static readonly ComplexNumber i = new ComplexNumber(0, 1);
        /// <summary>
        /// The component 'a' (Sometimes called 'x') is the real part of the complex number.<br></br>
        /// </summary>
        public float a;
        /// <summary>
        /// The component 'b' (Sometimes called 'y') is the imaginary part of the complex number.
        /// </summary>
        public float b;

        /// <summary>
        /// The component 'a' (Sometimes called 'x') is the real part of the complex number.
        /// The component 'b' (Sometimes called 'y') is the imaginary part of the complex number.
        /// </summary>
        public ComplexNumber(float a, float b)
        {
            this.a = a;
            this.b = b;
        }
        /// <summary>
        /// Makes a complex number with a distance from center of 1, and an angle of t.<br></br><br></br>
        /// Formula:<br></br>
        /// e^(it) = cos(t) + i * sin(t).
        /// </summary>
        public ComplexNumber(float t)
        {
            a = MathF.Cos(t);
            b = MathF.Sin(t);
        }
        #endregion

        #region Properties
        public static float Re(ComplexNumber z) => z.a;
        public static float Im(ComplexNumber z) => z.b;

        /// <summary>
        /// Returns a float representing the complex numbers distance from the origin.<br></br><br></br>
        /// Formula: <br></br>
        /// √a²+b².
        /// </summary>
        public float Magnitude() => MathF.Sqrt(a * a + b * b);
        public float SquareMagnitude() => a * a + b * b;
        /// <summary>
        /// Returns a float representing the angle of the complex number.<br></br><br></br>
        /// Formula: <br></br>
        /// arctan(y / x).
        /// </summary>
        public float Angle() => MathF.Atan2(b, a);

        /// <summary>
        /// Returns the distance between the origin and the complex number (r = √a²+b²).
        /// Also returns the angle of the complex number (theta = Atan2(b, a)).
        /// </summary>
        public (float r, float theta) ToPolar()
        {
            float r = Magnitude();
            float theta = Angle();

            return (r, theta);
        }
        #endregion

        #region Modifications
        /// <summary>
        /// Returns the conjugate of a complex number.<br></br>
        /// (Multiplies the imaginary component by (-1))
        /// </summary>
        public ComplexNumber Conj()
        {
            return new ComplexNumber(a, -b);
        }

        /// <summary>
        /// Swaps the values of a and b in a complex number.
        /// </summary>
        public ComplexNumber Swapped()
        {
            return new ComplexNumber(b, a);
        }

        /// <summary>
        /// Normalized a complex number, so its distance from the origin is 1.<br><br></br></br>
        /// Formula:<br></br>
        /// z / z.Magnitude().
        /// </summary>
        public ComplexNumber Normalize()
        {
            float mag = Magnitude();
            if (mag == 0)
                throw new DivideByZeroException();

            return new ComplexNumber(a, b) / mag;
        }
        /// <summary>
        /// Normalized a complex number, so its distance from the origin is 1.<br><br></br></br>
        /// Formula:<br></br>
        /// z / z.Magnitude().
        /// </summary>
        public ComplexNumber Copy()
        {
            return new ComplexNumber(a, b);
        }
        /// <summary>
        /// Normalized a complex number, so its distance from the origin is 1.<br><br></br></br>
        /// Formula:<br></br>
        /// z / z.Magnitude().
        /// </summary>
        public static void Normalize(ref ComplexNumber z)
        {
            float mag = z.Magnitude();
            if (mag == 0)
                throw new DivideByZeroException();

            z /= mag;
        }
        /// <summary>
        /// Normalized a complex number, so its distance from the origin is 1.<br><br></br></br>
        /// Formula:<br></br>
        /// z / z.Magnitude().
        /// </summary>
        public static ComplexNumber Normalize(ComplexNumber z)
        {
            float mag = z.Magnitude();
            if (mag == 0)
                throw new DivideByZeroException();

            return z / mag;
        }
        #endregion

        #region Operations

        #region Exponentiation
        /// <summary>
        /// Raises a complex number to any power.
        /// </summary>
        public ComplexNumber Pow(float power)
        {
            (float r, float theta) = ToPolar();

            float newA = MathF.Pow(r, power) * MathF.Cos(theta * power);
            float newB = MathF.Pow(r, power) * MathF.Sin(theta * power);
            return new ComplexNumber(newA, newB);
        }

        /// <summary>
        /// Raises a real number to a complex power.
        /// </summary>
        public static ComplexNumber Pow(float c, ComplexNumber complexPower)
        {
            float newA = MathF.Cos(complexPower.b * MathF.Log(c));
            float newB = MathF.Sin(complexPower.b * MathF.Log(c));
            return MathF.Pow(c, complexPower.a) * new ComplexNumber(newA, newB);
        }

        /// <summary>
        /// Raises a real number to a complex power.
        /// </summary>
        public static ComplexNumber Pow(int c, ComplexNumber complexPower)
        {
            float newA = MathF.Cos(complexPower.b * MathF.Log(c));
            float newB = MathF.Sin(complexPower.b * MathF.Log(c));
            return MathF.Pow(c, complexPower.a) * new ComplexNumber(newA, newB);
        }

        /// <summary>
        /// Raises a complex number to a complex power.<br></br><br></br>
        /// Formula:<br></br>
        /// (r₁ * e^(it₁))^(r₂ * e^(it₂)) = r₁^(r₂ * e^(it₂)) * e^(it₁ * r₂ * e^(it₂))<br></br>
        /// Multiplying e^(it₂) by i is the same as taking its conjugate and then swapping a and b.
        /// </summary>
        public static ComplexNumber Pow(ComplexNumber z, ComplexNumber complexPower)
        {
            (float r1, float t1) = z.ToPolar();
            (float r2, float t2) = complexPower.ToPolar();

            ComplexNumber rPart = Pow(r1, r2 * new ComplexNumber(t2));
            ComplexNumber ePart = Pow(MathF.E, t1 * r2 * new ComplexNumber(t2).Conj().Swapped());
            return rPart * ePart;
        }
        /// <summary>
        /// Raises a complex number to a complex power.<br></br><br></br>
        /// Formula:<br></br>
        /// (r₁ * e^(it₁))^(r₂ * e^(it₂)) = r₁^(r₂ * e^(it₂)) * e^(it₁ * r₂ * e^(it₂))
        /// </summary>
        public ComplexNumber Pow(ComplexNumber complexPower)
        {
            (float r1, float t1) = ToPolar();
            (float r2, float t2) = complexPower.ToPolar();

            ComplexNumber rPart = Pow(r1, r2 * new ComplexNumber(t2));
            // Multiplying e^(it₂) by i is the same as taking its conjugate and then swapping a and b.
            ComplexNumber ePart = Pow(MathF.E, t1 * r2 * new ComplexNumber(t2).Conj().Swapped());
            return rPart * ePart;
        }
        #endregion

        #region Trigonometry
        /// <summary>
        /// Returns the sine of a complex number.<br></br><br></br>
        /// Formula:<br></br>
        /// 1 / 2i * (e^(iz) - e^(-iz)).
        /// </summary>
        public static ComplexNumber Sin(ComplexNumber z)
        {
            return -0.5f * i * (Pow(MathF.E, z.Conj().Swapped()) - Pow(MathF.E, -z.Conj().Swapped()));
        }

        /// <summary>
        /// Returns the cosine of a complex number.<br></br><br></br>
        /// Formula:<br></br>
        /// 1 / 2 * (e^(i(z)) + e^(-i(z))).
        /// </summary>
        public static ComplexNumber Cos(ComplexNumber z)
        {
            return 0.5f * (Pow(MathF.E, z.Conj().Swapped()) + Pow(MathF.E, -z.Conj().Swapped()));
        }

        /// <summary>
        /// Returns the cosine of a complex number.<br></br><br></br>
        /// Formula:<br></br>
        /// sin(z) / cos(z).
        /// </summary>
        public static ComplexNumber Tan(ComplexNumber z)
        {
            return Sin(z) / Cos(z);
        }
        #endregion

        #region Logarithmic Functions
        /// <summary>
        /// Returns natural logarithm of a complex number.<br></br><br></br>
        /// Formula:<br></br>
        /// ln(r) + it.
        /// </summary>
        public static ComplexNumber Ln(ComplexNumber z)
        {
            float r = z.Magnitude();
            float theta = z.Angle();

            return new ComplexNumber(MathF.Log(r), theta);
        }

        /// <summary>
        /// Returns logarithm (with a real base) of a complex number.<br></br><br></br>
        /// Formula:<br></br>
        /// (ln(r) + it) / ln(base).
        /// </summary>
        public static ComplexNumber Log(float realBase, ComplexNumber z)
        {
            float r = z.Magnitude();
            float theta = z.Angle();

            return new ComplexNumber(MathF.Log(r), theta) / MathF.Log(realBase);
        }

        /// <summary>
        /// Returns logarithm (with a complex base) of a real number.<br></br><br></br>
        /// Formula:<br></br>
        /// ln(c) / ln(r) + it).
        /// </summary>
        public static ComplexNumber Log(ComplexNumber complexBase, float c)
        {
            float r = complexBase.Magnitude();
            float theta = complexBase.Angle();

            return MathF.Log(c) / new ComplexNumber(MathF.Log(r), theta);
        }

        /// <summary>
        /// Returns logarithm (with a complex base) of a complex number.<br></br><br></br>
        /// Formula:<br></br>
        /// (ln(r) + it) / ln(base).
        /// </summary>
        public static ComplexNumber Log(ComplexNumber complexBase, ComplexNumber z)
        {
            float r1 = complexBase.Magnitude();
            float theta1 = complexBase.Angle();
            float r2 = z.Magnitude();
            float theta2 = z.Angle();

            return new ComplexNumber(MathF.Log(r2), theta2) / new ComplexNumber(MathF.Log(r1), theta1);
        }

        /// <summary>
        /// Returns the product log of a complex number.<br></br><br></br>
        /// </summary>
        public static ComplexNumber ProductLog(ComplexNumber z)
        {
            ComplexNumber z_n = new(1, 1);
            ComplexNumber prevZ_n = new(-1, -1);

            Random rand = new Random();
            //Checks until the floats cant compute any more precisely
            while (Re(prevZ_n) != Re(z_n) && Im(prevZ_n) != Im(z_n))
            {
                prevZ_n = z_n.Copy();
                if (Re(z_n) != Re(z_n) || Im(z_n) != Im(z_n))
                {
                    //If x_n is NaN a new starting estimate is set (one higher than the previous)
                    //and this function will compute the product log and return it back to here
                    return ProductLog(z, new ComplexNumber((float)rand.NextDouble() * 2 * MathF.PI)
                        * (float)rand.NextDouble() * 0.25f);
                }
                z_n -= (z_n * Pow(MathF.E, z_n) - z) / ((z_n + 1) * Pow(MathF.E, z_n));
            }

            return z_n;
        }
        private static ComplexNumber ProductLog(ComplexNumber z, ComplexNumber startingEstimate)
        {
            ComplexNumber z_n = startingEstimate;
            ComplexNumber prevZ_n = -startingEstimate;
            if (startingEstimate.SquareMagnitude() == 0)
                startingEstimate = -new ComplexNumber(1, 1);

            Random rand = new Random();
            //Checks until the floats cant compute any more precisely
            while(Re(prevZ_n) != Re(z_n) && Im(prevZ_n) != Im(z_n))
            {
                prevZ_n = z_n.Copy();
                if (Re(z_n) != Re(z_n) || Im(z_n) != Im(z_n))
                {
                    //If x_n is NaN a new starting estimate is set (one higher than the previous)
                    //and this function will compute the product log and return it back to here
                    return ProductLog(z, new ComplexNumber((float)rand.NextDouble() * 2 * MathF.PI) 
                        * (float)rand.NextDouble() * (startingEstimate.Magnitude() + 0.25f));
                }
                z_n -= (z_n * Pow(MathF.E, z_n) - z) / ((z_n + 1) * Pow(MathF.E, z_n));
            }

            return z_n;
        }

        /// <summary>
        /// Returns the product log of a real number.<br></br><br></br>
        /// </summary>
        public static float ProductLog(float x)
        {
            float x_n = 1;
            float prevX_n = -1;

            //Checks until the floats cant compute any more precisely
            while (MathF.Abs(x_n - prevX_n) != 0) {
                if (x_n != x_n)
                {
                    //If x_n is NaN a new starting estimate is set (one higher than the previous)
                    //and this function will compute the product log and return it back to here
                    return ProductLog(x, 2);
                }
                prevX_n = x_n;
                x_n -=(x_n * MathF.Pow(MathF.E, x_n) - x) / (x_n * (x_n * MathF.Pow(MathF.E, x_n) + MathF.Pow(MathF.E, x_n)));
            }

            return x_n;
        }
        private static float ProductLog(float x, float startingEstimate)
        {
            float x_n = startingEstimate;
            float prevX_n = -startingEstimate;
            if (startingEstimate == 0)
                startingEstimate = -1;

            //Checks until the floats cant compute any more precisely
            while (MathF.Abs(x_n - prevX_n) != 0)
            {
                //If x_n is NaN the function is called recursively with a new starting estimate
                if (x_n != x_n)
                    return ProductLog(x, startingEstimate + 1);
                prevX_n = x_n;
                x_n -= (x_n * MathF.Pow(MathF.E, x_n) - x) / (x_n * (x_n * MathF.Pow(MathF.E, x_n) + MathF.Pow(MathF.E, x_n)));
            }

            return x_n;
        }
        #endregion

        #region Calculus

        /// <summary>
        /// Returns a very rough estimation of the factorial factorial of a complex number.<br></br><br></br>
        /// </summary>
        public static ComplexNumber Factorial(ComplexNumber input)
        {
            //double r = SpecialFunctions.Gamma(input);
            ComplexNumber[] d = new ComplexNumber[] { 
                new(2.48574089138753565546f * MathF.Pow(10, -5), 0),
                new(1.05142378581721974210f, 0),
                new(-3.45687097222016235469f, 0),
                new(4.51227709466894823700f, 0),
                new(-2.98285225323576655721f, 0),
                new(1.05639711577126713077f, 0),
                new(-1.95428773191645869583f * MathF.Pow(10, -1), 0),
                new(1.70970543404441224307f * MathF.Pow(10, -2), 0),
                new(-5.71926117404305781283f * MathF.Pow(10, -4), 0),
                new(4.63399473359905636708f * MathF.Pow(10, -6), 0),
                new(-2.71994908488607703910f * MathF.Pow(10, -9), 0)
            };
            ComplexNumber sum = d[0];
            for(int i = 1; i < 10; i++)
            {
                sum += d[i] / (input + i);
            }
            return 2 * MathF.Sqrt(MathF.E / MathF.PI) * Pow((input + 10.9005f + 0.5f) / MathF.E, input + 0.5f) * sum;
        }
        public static ComplexNumber ZetaFunction(ComplexNumber s, int iterations)
        {
            return EtaFunction(s, iterations) / (1 - Pow(2, 1-s));
        }
        public static ComplexNumber EtaFunction(ComplexNumber s, int iterations)
        {
            if (Re(s) <= 0)
                return new ComplexNumber(0, 0);

            ComplexNumber sum = new ComplexNumber(0, 0);
            for (int i = 1; i <= iterations; i++)
            {
                sum += MathF.Pow(-1, i + 1) * Pow(i, -s);
            }

            return sum;
        }

        //Defined for -∞ < Re(s) < 0 ∧ 1 < Re(s) < ∞
        public static ComplexNumber AnalyticZetaFunction(ComplexNumber s, int iterations)
        {
            if(Re(s) < 0)
                return Factorial(-s) * Pow(2 * MathF.PI, s - 1) * 2 * Cos(MathF.PI / 2 * (s - 1)) * ZetaFunction(1 - s, iterations);
            return ZetaFunction(s, iterations);
        }

        #endregion

        #region Complex Transformations
        /// <summary>
        /// Returns inverted value of a complex number.<br></br><br></br>
        /// </summary>
        public static ComplexNumber CircleInversion(ComplexNumber z)
        {
            return z / (z.a * z.a + z.b * z.b);
        }

        /// <summary>
        /// Returns inverted value of a complex number.<br></br><br></br>
        /// </summary>
        public static ComplexNumber CircleInversion(ComplexNumber z, ComplexNumber circleCenter)
        {
            ComplexNumber centerToZ = z - circleCenter;
            return centerToZ / (centerToZ.a * centerToZ.a + centerToZ.b * centerToZ.b);
        }

        /// <summary>
        /// Returns inverted value of a complex number.<br></br><br></br>
        /// </summary>
        public static ComplexNumber CircleInversion(ComplexNumber z, ComplexNumber circleCenter, float radius)
        {
            ComplexNumber centerToZ = z - circleCenter;
            return centerToZ * radius * radius / (centerToZ.a * centerToZ.a + centerToZ.b * centerToZ.b);
        }
        #endregion

        /// <summary>
        /// Returns a complex number as a string.
        /// Format: (a + bi)
        /// </summary>
        public override string ToString()
        {
            string outputString = "(";
            float roundedA = MathF.Round(a * 100000) / 100000;
            float roundedB = MathF.Round(b * 100000) / 100000;
            outputString += roundedA;

            if (roundedB == 0)
                return outputString + ")";

            if (roundedB >= 0)
                outputString += $" + {roundedB}i)";
            else
                outputString += $" - {-roundedB}i)";

            if (roundedA != 0)
                return outputString;


            if (roundedB >= 0)
                outputString = $"({roundedB}i)";
            else
                outputString = $"({-roundedB}i)";

            return outputString;
        }

        #region Standard Operations
        public static ComplexNumber operator +(ComplexNumber z, float c) => new ComplexNumber(z.a + c, z.b);
        public static ComplexNumber operator +(float c, ComplexNumber z) => new ComplexNumber(z.a + c, z.b);
        public static ComplexNumber operator +(ComplexNumber z, int c) => new ComplexNumber(z.a + c, z.b);
        public static ComplexNumber operator +(int c, ComplexNumber z) => new ComplexNumber(z.a + c, z.b);
        public static ComplexNumber operator +(ComplexNumber z1, ComplexNumber z2) => new ComplexNumber(z1.a + z2.a, z1.b + z2.b);
        public static ComplexNumber operator -(ComplexNumber z) => new ComplexNumber(-z.a, -z.b);
        public static ComplexNumber operator -(ComplexNumber z, float c) => new ComplexNumber(z.a - c, z.b);
        public static ComplexNumber operator -(float c, ComplexNumber z) => new ComplexNumber(c - z.a, -z.b);
        public static ComplexNumber operator -(ComplexNumber z, int c) => new ComplexNumber(z.a - c, z.b);
        public static ComplexNumber operator -(int c, ComplexNumber z) => new ComplexNumber(c - z.a, -z.b);
        public static ComplexNumber operator -(ComplexNumber z1, ComplexNumber z2) => new ComplexNumber(z1.a - z2.a, z1.b - z2.b);

        public static ComplexNumber operator *(ComplexNumber z, float c) => new ComplexNumber(z.a * c, z.b * c);
        public static ComplexNumber operator *(float c, ComplexNumber z) => new ComplexNumber(z.a * c, z.b * c);
        public static ComplexNumber operator *(ComplexNumber z, int c) => new ComplexNumber(z.a * c, z.b * c);
        public static ComplexNumber operator *(int c, ComplexNumber z) => new ComplexNumber(z.a * c, z.b * c);
        public static ComplexNumber operator *(ComplexNumber z1, ComplexNumber z2)
        {
            (float r1, float theta1) = z1.ToPolar();
            (float r2, float theta2) = z2.ToPolar();

            float a = r1 * r2 * MathF.Cos(theta1 + theta2);
            float b = r1 * r2 * MathF.Sin(theta1 + theta2);

            return new ComplexNumber(a, b);
        }

        public static ComplexNumber operator /(ComplexNumber z, float c) => new ComplexNumber(z.a / c, z.b / c);
        public static ComplexNumber operator /(float c, ComplexNumber z) => c * z.Pow(-1);
        public static ComplexNumber operator /(ComplexNumber z, int c) => new ComplexNumber(z.a / c, z.b / c);
        public static ComplexNumber operator /(int c, ComplexNumber z) => c * z.Pow(-1);
        public static ComplexNumber operator /(ComplexNumber z1, ComplexNumber z2)
        {
            (float r1, float theta1) = z1.ToPolar();
            (float r2, float theta2) = z2.ToPolar();

            float a = r1 / r2 * MathF.Cos(theta1 - theta2);
            float b = r1 / r2 * MathF.Sin(theta1 - theta2);

            return new ComplexNumber(a, b);
        }

        public static ComplexNumber operator %(ComplexNumber z, float c) => new ComplexNumber(z.a % c, z.b % c);
        public static ComplexNumber operator %(float c, ComplexNumber z) => new ComplexNumber(z.a % c, z.b % c);
        public static ComplexNumber operator %(ComplexNumber z, int c) => new ComplexNumber(z.a % c, z.b % c);
        public static ComplexNumber operator %(int c, ComplexNumber z) => new ComplexNumber(z.a % c, z.b % c);
        public static ComplexNumber operator %(ComplexNumber z1, ComplexNumber z2) => new ComplexNumber(z1.a % z2.a, z1.b % z2.b);

        public static ComplexNumber operator ++(ComplexNumber z1) => new ComplexNumber(z1.a++, z1.b++);
        public static ComplexNumber operator --(ComplexNumber z1) => new ComplexNumber(z1.a--, z1.b--);

        public static explicit operator ComplexNumber((float r, float theta)polar) 
            => new ComplexNumber(polar.r * MathF.Cos(polar.theta), polar.r * MathF.Sin(polar.theta));
        #endregion
        #endregion
    }
}

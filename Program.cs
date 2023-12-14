using System;
using System.IO;
using System.Reflection;
using System.Text;
using System.Threading;

namespace ComplexNumbers
{
    class Program
    {
        private static ComplexNumber XiFunction(float t)
        {
            ComplexNumber s = new ComplexNumber(0.5f, t);
            ComplexNumber zeta = ComplexNumber.ZetaFunction(s, 2000);
            //zeta *= ComplexNumber.Pow(MathF.PI, -s / 2f);
            //zeta *= (s - 1);
            //zeta *= ComplexNumber.Factorial(s / 2f);
            return zeta;
        }

        static void Main(string[] args)
        {
            //Xi function:
            Console.WriteLine($"Xi: {XiFunction(15f)}");

            Thread searchThread = new Thread(SearchThread);
            searchThread.Start();

        }

        static float epsilon = 0.1f;
        static void SearchThread() //Find potential zeros of the zeta function
        {
            ComplexNumber previousValue = XiFunction(0f);
            int index = 0;
            float stepSize = 0.002f;
            int roots = 0;
            StringBuilder stringBuilder = new StringBuilder();

            while (true){
                index++;
                ComplexNumber currentValue = XiFunction(index * stepSize);
                if (previousValue.Magnitude() <= epsilon &&
                    MathF.Sign(ComplexNumber.Im(previousValue))
                    != MathF.Sign(ComplexNumber.Im(currentValue)))
                {
                    roots++;
                    Console.WriteLine($"Root #{roots} between: xi({(index - 1) * stepSize}) " +
                        $"and xi({index * stepSize})");
                    Console.WriteLine($"zeta(1/2 + {(index - 1) * stepSize}i) = {previousValue}");
                    Console.WriteLine($"zeta(1/2 + {index * stepSize}i) = {currentValue}");
                    Console.WriteLine();
                    stringBuilder.Append($"Root #{roots} between: xi({(index - 1) * stepSize}) " +
                        $"and xi({index * stepSize})\n");
                    File.WriteAllText(@"C:\Users\noah0\source\repos\OpenGL Math\zeta-roots\Zeros.txt", 
                        stringBuilder.ToString());
                    Thread.Sleep(0);
                }
                previousValue = currentValue.Copy();
            }
        }
    }
}

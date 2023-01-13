#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double f1(double x)
{
  //#  return /* Hier f_1(x) */;
  // x^3 -2x+2
  return pow(x, 3) - 2 * x + 2;
}

double df1(double x)
{
  // 3x^2-2
  return 3 * pow(x, 2) - 2;
}

double f2(double x)
{
  // sin((x*pi)/2) - cos((x*pi)*e^-3x
  return sin((x * M_PI) / 2) - cos(x * M_PI) * exp(-3 * x);
}

double df2(double x)
{
  double h = pow(10, -3);

  return ((f2(x + h) - f2(x - h)) / (2 * h));
}

double calculateBisection(double (*f)(double), double l, double r)
{
  static int n = 1;
  double fl = f(l);
  double fr = f(r);
  printf("  %.02d: %.12lf   %.12lf   %.12lf %.12lf\n", n, l, r, fl, fr);

  double m = 0.5 * (l + r);
  double fm = f(m);

  n++;
  if (fabs(r - l) <= pow(10, -12))
  {
    n = 1;
    return m;
  }
  else if (fm * fr > 0)
  {
    return calculateBisection(f, l, m);
  }
  else
  {
    return calculateBisection(f, m, r);
  }
}

double Bisektion(double (*f)(double), double xa, double xb)
{
  // Kopfzeile für die Ausgabe
  printf("  n:  l                 r                 f(l)               f(r)\n");

  return calculateBisection(f, xa, xb);
}

double Sekanten(double (*f)(double), double x0, double x1)
{
  static int n = 1;
  double fx0 = f(x0);
  double fx1 = f(x1);
  double x2 = x0 - fx0 * ((x1 - x0) / (fx1 - fx0));

  printf("%.02d: x0: %.12lf, x1: %.12lf, f(x0): %.12lf, f(x1): %12.lf, x2: %.12lf \n", n, x0, x1, fx0, fx1, x2);
  if (fabs(x2 - x1) < pow(10, -12))
  {
    n = 1;
    return x2;
  }

  n++;
  return Sekanten(f, x1, x2);
}

double calculateNewton(double (*f)(double), double (*dfdx)(double), double x0)
{
  static int n = 1;
  double fx0 = f(x0);
  double dfdxX0 = dfdx(x0);
  double x1 = x0 - (fx0 / dfdxX0);

  printf("%.03d: %16.12lf %16.12lf %16.12lf %16.12lf \n", n, x0, fx0, dfdxX0, x1);
  n++;

  if (fabs(x1 - x0) < pow(10, -12))
  {
    n = 1;
    return x1;
  }
  else if (n > 100)
  {
    n = 1;
    return NAN;
  }

  return calculateNewton(f, dfdx, x1);
}

double Newton(double (*f)(double), double (*dfdx)(double), double x0)
{
  printf("  n: x0                f(x0)            f'(x0)           x1\n");
  return calculateNewton(f, dfdx, x0);
}

int main()
{
  double result;

  printf("\n\nTask 1:\n");

  printf("\n\nBisektion Start -3.0, -0.5\n–––––––––\n");
  result = Bisektion(f1, -3.0, -0.5);
  printf("Nullstelle Bisektion: %.14lf\n", result);

  printf("\n\nSekanten Start -3.0, -0.5\n––––––––\n");
  result = Sekanten(f1, -3.0, -0.5);
  printf("Nullstelle Sekanten: %.14lf\n", result);

  printf("\n\nNewton Start -0.5\n––––––\n");
  result = Newton(f1, df1, -0.5);
  printf("Nullstelle: %.14lf\n", result);

  printf("\n\nNewton Start -0.51\n––––––\n");
  result = Newton(f1, df1, -0.51);
  printf("Nullstelle: %.14lf\n", result);

  printf("\n\nTask 2:\n");
  printf("\n\nBisektion Start 0.0, 3.0\n–––––––––\n");
  result = Bisektion(f2, 0.0, 3.0);
  printf("Nullstelle Bisektion: %.14lf\n", result);

  printf("\n\nSekanten Start 0.0, 1.0\n––––––––\n");
  result = Sekanten(f2, 0.0, 1.0);
  printf("Nullstelle Sekanten: %.14lf\n", result);

  printf("\n\nNewton Start 0.9\n––––––\n");
  result = Newton(f2, df2, 0.9);
  printf("Nullstelle: %.14lf\n", result);

  return 0;
}

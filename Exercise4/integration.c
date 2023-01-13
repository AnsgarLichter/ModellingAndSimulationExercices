#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/**
 * @brief fa(x)=sic(x)
 *
 * @param x x
 *
 * @return f(x)
 */

double f1(double x)
{
  if (x == 0)
  {
    return 1;
  }
  return sin(x) / x;
}

double f2(double x)
{
  return exp(x);
}

/**
 * @brief Integration von f mit Untersummen.
 *
 * @param a untere Grenze des Integrationsbereiches
 * @param b obere Grenze des Integrationsbereiches
 * @param n Anzahl der Stützstellen
 * @param f zu integrierende Funktion
 */

double int_unter(double a, double b, long n, double (*f)(double))
{
  double h = (b - a) / n;
  double sum = 0;
  for (int i = 0; i < n; i++)
  {
    sum += f(a + i * h);
  }
  return h * sum;
}

/**
 * @brief Integration von f mit Mittelsummen.
 *
 * @param a untere Grenze des Integrationsbereiches
 * @param b obere Grenze des Integrationsbereiches
 * @param n Anzahl der Stützstellen
 * @param f zu integrierende Funktion
 */

double int_mitte(double a, double b, long n, double (*f)(double))
{
  double h = (b - a) / n;
  double intervalmitte = 0.0;
  double sum = 0.0;

  for (int i = 0; i < n; i++)
  {
    intervalmitte = 0.5 * ((a + (i + 1) * h) + (a + i * h));
    sum += f(intervalmitte);
  }

  return h * sum;
}

/**
 * @brief Integration von f mit Obersummen.
 *
 * @param a untere Grenze des Integrationsbereiches
 * @param b obere Grenze des Integrationsbereiches
 * @param n Anzahl der Stützstellen
 * @param f zu integrierende Funktion
 */

double int_ober(double a, double b, long n, double (*f)(double))
{
  double h = (b - a) / n;
  double sum = 0.0;

  for (int i = 0; i < n; i++)
  {
    sum += f(a + (i + 1) * h);
  }
  return h * sum;
}

/**
 * @brief Integration von f mit der Trapez-Regel.
 *
 * @param a untere Grenze des Integrationsbereiches
 * @param b obere Grenze des Integrationsbereiches
 * @param n Anzahl der Stützstellen
 * @param f zu integrierende Funktion
 */

double int_trapez(double a, double b, long n, double (*f)(double))
{
  return 0.5 * (int_unter(a, b, n, f) + int_ober(a, b, n, f));
}

/**
 * @brief Integration von f mit der Simpson-Regel.
 *
 * @param a untere Grenze des Integrationsbereiches
 * @param b obere Grenze des Integrationsbereiches
 * @param n Anzahl der Stützstellen
 * @param f zu integrierende Funktion
 */

double int_simpson(double a, double b, long n, double (*f)(double))
{
  double h = (b - a) / n;
  double sum = 0.0;

  for (int i = 0; i < n; i++)
  {
    double xi = a + i * h;
    double xi1 = a + (i + 0.5) * h;
    double xi2 = a + (i + 1) * h;

    sum += (h / 6.0) * (f(xi) + 4 * f(xi1) + f(xi2));
  }

  return sum;
}

/**
 * @brief Integration von f mit der Trapez-Regel.
 *
 * @param a untere Grenze des Integrationsbereiches
 * @param b obere Grenze des Integrationsbereiches
 * @param n Anzahl der Stützstellen
 * @param f zu integrierende Funktion
 */

double bogenlaenge(double a, double b, long n, double (*f)(double))
{
  double h = ((b - a) / n);
  double length = 0;

  for (int i = 0; i < n; i++)
  {
    double xi = a + i * h;
    double xi1 = a + (i + 1) * h;
    length += sqrt(pow(h, 2) + pow((f(xi1) - f(xi)), 2));
  }

  return length;
}

/**
 * @brief Sucht die Anzahl der Stützstellen, die benötigt werden um einen Approximationsfehler < err zu haben.
 *
 * @param a        untere Grenze des Integrationsbereiches
 * @param b        obere Grenze des Integrationsbereiches
 * @param f        zu integrierende Funktion
 * @param integral zu verwendendes Integrationsverfahren
 * @param real     analytische Lösung
 * @param err      Fehlertoleranz
 */

void findsteps(double a, double b, double (*f)(double), double integral(double, double, long, double (*)(double)), double real, double err)
{
  int n = 0;
  double error;

  do
  {
    n++;
    error = real - integral(a, b, n, f);
  } while (err <= fabs(error));

  printf("%6d Stützstellen, err=%le\n", n, error);
}

int main(void)
{
  int n;
  double real;
  double a = 0.0;
  double b = 2.0;
  double (*f)(double);

  // Sie können f benutzen und an dieser Stelle abändern wenn Sie fb haben wollen
  f = f1;
  real = 1.6054129768026948485767201; // b=2

  // f = f2;
  //  real = exp(b)- exp(a);

  for (n = 4; n <= 128; n = n * 2)
  {
    double unten = int_unter(a, b, n, f);
    double mitte = int_mitte(a, b, n, f);
    double oben = int_ober(a, b, n, f);
    double trapez = int_trapez(a, b, n, f);
    double simpson = int_simpson(a, b, n, f);

    printf("Stützstellen   h        Simpson  Trapez   Oben     Mitte    Unten\n");
    printf("%.03d            %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf\n", n, ((b - a) / n), simpson, trapez, oben, mitte, unten);
    printf("%.03d            %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf\n\n\n", n, ((b - a) / n), real - simpson, real - trapez, real - oben, real - mitte, real - unten);
  }

  double err = 1e-3;

  printf("\n Rechteck Untere Summe: \n");
  findsteps(a, b, f, int_unter, real, err);
  printf("\n Rechteck Mittlere Summe: \n");
  findsteps(a, b, f, int_mitte, real, err);
  printf("\n Rechteck Obere Summe: \n");
  findsteps(a, b, f, int_ober, real, err);
  printf("\n Trapez Summe: \n");
  findsteps(a, b, f, int_trapez, real, err);
  printf("\n Simpson Summe: \n");
  findsteps(a, b, f, int_simpson, real, err);

  printf("\n\n");
  for (int i = 3; i < 9; i++)
  {
    err = pow(10, -i);

    printf("\ne^-%02d:\n", i);

    printf("\n Rechteck Mittlere Summe: \n");
    findsteps(a, b, f, int_mitte, real, err);
    printf("\n Rechteck Obere Summe: \n");
    findsteps(a, b, f, int_ober, real, err);
    printf("\n Simpson Summe: \n");
    findsteps(a, b, f, int_simpson, real, err);
  }

  printf("\n Bogenlänge: \n");
  double bogenlaengeValue = bogenlaenge(a, b, n, f2);

  printf(" %.6lf \n", bogenlaengeValue);
  return 0;
}

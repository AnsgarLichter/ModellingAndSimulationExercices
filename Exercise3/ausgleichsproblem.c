#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*********************************************************************/
/*                                                                   */
/* Es liegen Messdaten "input.dat" in Form von Wertepaaren vor.      */
/* Durch diese Punkte muss eine Ausgleichskurve gelegt werden.       */
/*                                                                   */
/*********************************************************************/
/*             Zu Beginn einige Hilfsfunktionen                      */
// --------------------------------------------------------------------

/*
------- AT*A ---------
2.500000e+01 3.000000e+02
3.000000e+02 4.900000e+03
------- AT*logI -----
8.282651e+01
1.223997e+03
------- ATAinv ------
1.507692e-01 -9.230769e-03
-9.230769e-03 7.692308e-04
------- lambda -----
ln_a: 1.189257e+00
b:    1.769836e-01
*/

/** @brief Funktion um die Wertepaare abzuzählen.
 */
int getNumberofPoints(char *name)
{
  FILE *fp;
  char *line = NULL;
  size_t len = 0;

  if ((fp = fopen(name, "r")) == NULL)
  {
    exit(EXIT_FAILURE);
  }

  int cnt = 0;
  while (getline(&line, &len, fp) != -1)
  {
    cnt++;
  }

  free(line);
  fclose(fp);

  return cnt;
}
/** @brief In dieser Funktion werden die Wertepaare eingelesen und in Form von Arrays t[N] und I[N] übergeben.
 */
void readFile(char *name, int t[], double I[])
{
  FILE *fp;
  char *line = NULL;
  size_t len = 0;

  if ((fp = fopen(name, "r")) == NULL)
  {
    exit(EXIT_FAILURE);
  }

  int cnt = 0;
  while (getline(&line, &len, fp) != -1)
  {
    sscanf(line, "%i %lf", &t[cnt], &I[cnt]);
    cnt++;
  }

  free(line);
  fclose(fp);
}

// --------------------------------------------------------------------
// --------------------------------------------------------------------

int main(int argc, char *argv[])
{
  // Abzaehlen der Wertepaare
  int N = getNumberofPoints("input.dat");

  int t[N];    // Vektor für den Kalenderwoche der Messung
  double I[N]; // Vektor für den gemessenen Inzidenzwerte
  // Einlesen der Daten
  readFile("input.dat", t, I);

  // Hilsvariabelen fuer t_0 und Logarithmen von I
  int t_0 = t[0];

  double logI[N];

  // Hilfsvektor y erstellen:
  printf("N: %03d \n", N);
  for (int i = 0; i < N; i++)
  {
    logI[i] = log(I[i]);
  }

  // A & AT erstellen
  double A[N][2];
  double AT[2][N];

  for (int i = 0; i < N; i++)
  {
    double tDifference = t[i] - t_0;

    A[i][0] = 1;
    A[i][1] = tDifference;

    AT[0][i] = 1;
    AT[1][i] = tDifference;
  }

  printf("\nA:\n");
  for (int i = 0; i < N; ++i)
  {
    for (int j = 0; j < 2; ++j)
    {
      printf("%2.0lf  ", A[i][j]);
    }
    printf("\n");
  }

  printf("\nAT:\n");
  for (int i = 0; i < 2; ++i)
  {
    for (int j = 0; j < N; ++j)
    {
      printf("%2.0lf  ", AT[i][j]);
    }
    printf("\n");
  }

  //_____________________________________________________________________
  // Berechnung der Matrix A^T*A und der rechten Seiten A^T*y
  double ATA[2][2] = {{0.0, 0.0},
                      {0.0, 0.0}};
  double ATlogI[2] = {0.0, 0.0};

  // ATA berechnen
  printf("\nCalculate ATA\n");
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      for (int k = 0; k < N; k++)
      {
        ATA[i][j] += AT[i][k] * A[k][j];
      }
    }
  }

  printf("\nATA:\n");
  for (int i = 0; i < 2; ++i)
  {
    for (int j = 0; j < 2; ++j)
    {
      printf("%10.6le  ", ATA[i][j]);
    }
    printf("\n");
  }

  // AT*y berechnen
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < N; j++)
    {
      ATlogI[i] += AT[i][j] * logI[j];
    }
  }

  printf("\nATlogI:\n");
  for (int i = 0; i < 2; ++i)
  {
    printf("%10.6le\n", ATlogI[i]);
  }

  // Invertierung von ATA (wegen kleiner Groesse)
  double ATAInverted[2][2] = {{0.0, 0.0},
                              {0.0, 0.0}};

  double inversion_factor = 1.0 / ((ATA[0][0] * ATA[1][1]) - (ATA[0][1] * ATA[1][0]));
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      if (i != j)
      {
        ATAInverted[i][j] = inversion_factor * (-1) * ATA[i][j];
      }
      else
      {
        ATAInverted[1 - i][1 - j] = inversion_factor * ATA[i][j];
      }
    }
  }

  printf("\nATA Inverted:\n");
  for (int i = 0; i < 2; ++i)
  {
    for (int j = 0; j < 2; ++j)
    {
      printf("%10.6le  ", ATAInverted[i][j]);
    }
    printf("\n");
  }

  // Ermittlung von ln(a) und b aus ATA^{-1} und ATy
  double ln_a = 0;
  double b = 0;

  for (int i = 0; i < 2; i++)
  {
    ln_a += ATAInverted[0][i] * ATlogI[i];
    b += ATAInverted[1][i] * ATlogI[i];
  }

  printf("\nResult for a and b:\n");
  printf("ln_a = %.6le \nb = %.6le\n", ln_a, b);

  // Bestimmung von a selbst
  double a = exp(ln_a);

  // Plotten wenn plotflag!=0
  long plotflag = 1;

  if (plotflag)
  {
    FILE *gp = popen("gnuplot -p", "w");
    fprintf(gp, "reset; set key left top box; set xlabel \"t - t_0\";\n"
                "set ylabel \"y\";\n"
                "set autoscale fix\n"
                "set logscale xy\n"
                " f(t) = %le*exp(%le*(t - 23));\n "
                "set terminal png size 1920,1080; set output 'xyz.png'; plot f(x) lt -1 lw 2, \"input.dat\" using 1:2 pt 7 title 'measured data';\n", // lt: LineType, lw: Linewidth, using 1:2: zweite spalte verwenden
            a, b);
    pclose(gp);
  }

  return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define T_OFEN 6.0
#define WAERMEFLUSS 2.0
#define WAERMEQUELLE 0.6
#define L 5

#define time 100
#define DELTA_X 0.5
#define DELTA_T 0.1

/*
 *
-2.000000e+00 2.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
-4.000000e+00 8.000000e+00 -4.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
0.000000e+00 -4.000000e+00 8.000000e+00 -4.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
0.000000e+00 0.000000e+00 -4.000000e+00 8.000000e+00 -4.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
0.000000e+00 0.000000e+00 0.000000e+00 -4.000000e+00 8.000000e+00 -4.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 -4.000000e+00 8.000000e+00 -4.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 -4.000000e+00 8.000000e+00 -4.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 -4.000000e+00 8.000000e+00 -4.000000e+00 0.000000e+00 0.000000e+00
0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 -4.000000e+00 8.000000e+00 -4.000000e+00 0.000000e+00
0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 -4.000000e+00 8.000000e+00 -4.000000e+00
0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00
Fehler erste variante |T_1 - T_exakt|: 9.706634e-01

8.000000e+00 -8.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
-4.000000e+00 8.000000e+00 -4.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
0.000000e+00 -4.000000e+00 8.000000e+00 -4.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
0.000000e+00 0.000000e+00 -4.000000e+00 8.000000e+00 -4.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
0.000000e+00 0.000000e+00 0.000000e+00 -4.000000e+00 8.000000e+00 -4.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 -4.000000e+00 8.000000e+00 -4.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 -4.000000e+00 8.000000e+00 -4.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 -4.000000e+00 8.000000e+00 -4.000000e+00 0.000000e+00 0.000000e+00
0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 -4.000000e+00 8.000000e+00 -4.000000e+00 0.000000e+00
0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 -4.000000e+00 8.000000e+00 -4.000000e+00
0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00
Fehler zweite variante |T_2 - T_exakt|: 1.087792e-15
 *
 */

/**
 * @brief Temperaturfeld in Datei schreiben
 *
 * @param T        Temperaturfeld
 * @param timestep Zeitschritt
 * @param file     Ausgabedatei
 */
void writeDatFile(double *T, int timestep, FILE *file)
{
  int NX = (int)L / DELTA_X + 1;
  int k;

  for (k = 0; k < NX; k++)
  {
    fprintf(file, "%f %d %lf\n", k * DELTA_X, timestep, T[k]);
  }
  fprintf(file, "\n");
}

void printTridiagonalMatrix(double *b, double *d, double *c, int N)
{
  for (int i = 0; i < N; i++)
  {
    if (i > 0)
    {
      for (int j = 0; j < i - 1; j++)
      {
        printf("%le ", 0.0);
      }
      printf("%le ", b[i - 1]);
    }
    printf("%le ", d[i]);
    if (i < N - 1)
    {
      printf("%le ", c[i]);
      for (int j = i + 2; j < N; j++)
      {
        printf("%le ", 0.0);
      }
    }
    printf("\n");
  }
  printf("\n");
}

/**
 * @brief Loesung eines tridigaonalen Gleichungssystems
 *
 * @param b     Vektor mit Eintraegen der unteren Nebendiagonalen
 * @param d     Vektor mit Einträgen der Diagonalen
 * @param c     Vektor mit Einträgen der oberen Nebendiagonalen
 * @param x     Loesungsvektor
 * @param p     rechte Seite
 * @param N     Anzahl der Gleichungen
 */
void solveTridiagonalSystem(double *b, double *d, double *c, double *x, double *p, int N)
{
  double beta;
  double c_bar[N];

  // Vorwärts-Elimination
  x[0] = p[0];
  beta = d[0];
  for (int i = 0; i < N - 1; i++)
  {
    c_bar[i] = c[i] / beta;
    x[i] /= beta;
    beta = d[i + 1] - b[i] * c_bar[i];
    x[i + 1] = p[i + 1] - b[i] * x[i];
  }
  x[N - 1] /= beta;

  // Rückwärtssubstitution
  for (int i = N - 2; i >= 0; i--)
  {
    x[i] -= c_bar[i] * x[i + 1];
  }
}

double calcT_comp(double x)
{
  return T_OFEN - WAERMEFLUSS * (L - x) + 0.5 * WAERMEQUELLE * (L * L - x * x);
}

double calcNormDiffT(double *T, int N)
{
  double T_comp, sum;

  T_comp = calcT_comp(0.0);
  sum = 0.5 * (T[0] - T_comp) * (T[0] - T_comp);
  for (int i = 1; i < N - 1; i++)
  {
    T_comp = calcT_comp(i * DELTA_X);
    sum += (T[i] - T_comp) * (T[i] - T_comp);
  }
  T_comp = calcT_comp(L);
  sum += 0.5 * (T[N - 1] - T_comp) * (T[N - 1] - T_comp);

  return sqrt(sum * DELTA_X);
}

/**
 * @brief Berechnung der Wärmeleitung mit dem expliziten Euler-Verfahren
 *
 * @param Tneu     neues Temperaturfeld
 * @param T        aktuelles Temperaturfeld
 */
void calc(double *Tneu, double *T)
{
  int NX = (int)L / DELTA_X + 1;
  Tneu[NX - 1] = T_OFEN;

  double deltaXSquare = DELTA_X * DELTA_X;
  for (int i = 1; i <= NX - 2; i++)
  {
    Tneu[i] = T[i] + DELTA_T * (WAERMEQUELLE + (T[i + 1] - 2 * T[i] + T[i - 1]) / (deltaXSquare));
  }

  Tneu[0] = T[1] - WAERMEFLUSS * DELTA_X;
}

int main()
{
  int NX = (int)L / DELTA_X + 1;
  int TIMESTEPS = (int)time / DELTA_T;
  // Ausgabedatei
  FILE *heatfile;
  heatfile = fopen("heat.dat", "w");

  // Temperaturfeld
  double *T = (double *)calloc(NX, sizeof(double));

  // Konstruktion Gleichungssystem
  double *f = (double *)calloc(NX, sizeof(double));
  double *d = (double *)calloc(NX, sizeof(double));     // Diagonale
  double *b = (double *)calloc(NX - 1, sizeof(double)); // untere Nebendiagonale
  double *c = (double *)calloc(NX - 1, sizeof(double)); // obere Nebendiagonale

  double deltaXPow = pow(DELTA_X, 2);
  for (int i = 1; i < NX - 1; i++)
  {
    f[i] = WAERMEQUELLE;
    c[i] = -1.0 / deltaXPow;
    d[i] = 2.0 / deltaXPow;
  }

  for (int i = 0; i < NX - 2; i++)
  {
    b[i] = -1.0 / deltaXPow;
  }

  // Randbedingungen am rechten Rand
  b[NX - 2] = 0;
  d[NX - 1] = 1;
  f[NX - 1] = T_OFEN;

  // 1. Variante Randbedingung links
  f[0] = WAERMEFLUSS;
  c[0] = 1.0 / DELTA_X;
  d[0] = -1.0 / DELTA_X;

  printTridiagonalMatrix(b, d, c, NX);
  solveTridiagonalSystem(b, d, c, T, f, NX);
  printf("Fehler erste variante |T_1 - T_exakt|: %e\n", calcNormDiffT(T, NX));

  // 2. Variante Randbedingung links
  c[0] = -1.0 / (deltaXPow / 2.0);
  d[0] = 1.0 / (deltaXPow / 2.0);
  f[0] = WAERMEQUELLE - (WAERMEFLUSS / (DELTA_X / 2.0));

  printTridiagonalMatrix(b, d, c, NX);
  solveTridiagonalSystem(b, d, c, T, f, NX);
  printf("Fehler zweite variante |T_2 - T_exakt|: %e\n", calcNormDiffT(T, NX));

  free(f);
  free(d);
  free(b);
  free(c);

  // Berechnung der zeitlichen Entwicklung

  // Speicher für neue Temperatur
  double *Tneu = (double *)calloc(NX, sizeof(double));

  // Anfangsbedingung
  T[NX - 1] = T_OFEN;

  for (int i = 1; i < NX - 1; i++)
  {
    T[i] = WAERMEQUELLE;
  }

  T[0] = T[1] - WAERMEFLUSS * DELTA_X;



  for (int t = 0; t * DELTA_T < time + 1; t++)
  {
    if (t % (TIMESTEPS / 50) == 0)
    {
      writeDatFile(T, t * DELTA_T, heatfile);
    }

    // Berechnung des Feldes im nächsten Zeitschritt (Tneu)
    calc(Tneu, T);

    // Neu alt tauschen
    double *temp = Tneu;
    Tneu = T;
    T = temp;
  }

  // Aufräumen
  free(T);
  free(Tneu);
  fclose(heatfile);

  long plotflag = 1;

  if (plotflag)
  {
    char plotbefehl[1000];
    sprintf(plotbefehl, "echo "
                        "'set pm3d;\n"
                        "set palette defined (0 \"blue\", 1 \"white\", 4 \"yellow\", 6 \"red\");\n"
                        "set xlabel \"Raum\";\n"
                        "set ylabel \"Zeit\";\n"
                        "set zlabel \"Temperatur\";\n"
                        "splot \"heat.dat\"  every 1 w l lt 7 notitle;\n"
                        "' | gnuplot -persist");
    system(plotbefehl);
  }
  return 0;
}

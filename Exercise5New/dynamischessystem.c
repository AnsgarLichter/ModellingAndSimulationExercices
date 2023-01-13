#include "dynamischessystem.h"

/** Die rechten Seiten der Differentialgleichungen werden unter dem Pointer @*dgl abgelegt */
typedef double (*dgl)(double, double[2]);

/**
 * @brief Stellt die rechte Seite der Lotka-Volterra-Gleichung zur Berechnung des Räuber-Verhaltens dar.
 *
 * @param t   Zeitpunkt der Auswertung
 * @param x   Vector bestehend aus Zahl der Räuber- und Beutetiere
 */
double raeuber(double t, double x[2])
{
  return (-1.1 + 0.02 * x[1]) * x[0];
}

/**
 * @brief Stellt die rechte Seite der Lotka-Volterra-Gleichung zur Berechnung des Beute-Verhaltens dar.
 *
 * @param t   Zeitpunkt der Auswertung
 * @param x   Vector bestehend aus Zahl der Räuber- und Beutetiere
 */
double beute(double t, double x[2])
{
  return (0.35 - 0.03 * x[0]) * x[1];
}

/** Im Folgenden wird die Lotka-Volterra-Gleichungen numerisch integriert,
 *  dazu wurde bereits Runge-Kutta-Verfahren 2. Ordnung implementiert.
 *  Sie müssen noch folgende Verfahren implementieren
 *
 *    1.) Euler Verfahren
 *    2.) Runge-Kutta-Verfahren 4. Ordnung
 *
 *  Alle Verfahren haben folgende Argumente
 *
 * @param t     Zeit
 * @param dt    Zeitschrittweite
 * @param x     Vektor mit Räuber- und Beutetieren
 * @param fa    die 1. DGL
 * @param fb    die 2. DGL
 *
 * und sind unter dem Zeiger method abgelegt:
 */
typedef void (*method)(double, double, double[2], dgl, dgl);

/**  Die Verfahren werden auf die Lotka-Volterra-Gleichungen angewendet,
 *   die zeitliche Iteration wird mit der vorimplentierten Funktion @solve ausgeführt.
 *
 * @param dt        Zeitschrittweite
 * @param t_ende    Zeitende
 * @param x         Anfangswerte für die Räuber- und Beutetiere
 * @param raeuber   Räuber DGL
 * @param beute     Beute DGL
 * @param verfahren benutzte Integrationsmethode
 * @param file      Dateiname für die Ausgabe
 */
void solve(double dt, double t_ende, double x_start[2], dgl raeuber, dgl beute, method verfahren, char *file)
{
  FILE *plot = fopen(file, "w");
  double t = 0;
  // kopieren der Anfangswerte um Überschreiben zu vermeiden
  double x[2];
  x[0] = x_start[0];
  x[1] = x_start[1];

  // Rausschreiben des Anfangswerts
  plotCoordinate(plot, x[0], x[1], t);

  while (t < t_ende)
  {
    t += dt;
    verfahren(t, dt, x, raeuber, beute);

    // Rausschreiben jedes 100.ten Iterationsschritt
    if ((int)(ceil(t / dt)) % 100 == 0)
    {
      plotCoordinate(plot, x[0], x[1], t);
    }
  }
  fclose(plot);
  addPlot(gp, file);
}

/**
 * @brief Stellt einen Schritt des Runge-Kutta-Verfahrens zweiter Ordnung zur Lösung
 *        des Anfangswertproblems eines dynamischen Systems mit zwei Gleichungen dar.
 *
 * @param t     Zeit
 * @param dt    Zeitschrittweite
 * @param x     Zustandsvektor zum aktuellen Zeitpunkt
 * @param RT    die 1. DGL
 * @param BT    die 2. DGL
 *
 * @return  Zustandsvektor zum neuen Zeitschritt
 */
void runge_kutta2(double t, double dt, double x[2], dgl RT, dgl BT)
{
  double xtmp[2];

  xtmp[0] = x[0] + 0.5 * dt * RT(t, x);
  xtmp[1] = x[1] + 0.5 * dt * BT(t, x);

  x[0] += dt * RT(t + 0.5 * dt, xtmp);
  x[1] += dt * BT(t + 0.5 * dt, xtmp);
}

/**
 * @brief  Stellt einen Schritt des Euler-Verfahrens zur Lösung
 *         des Anfangswertproblems eines dynamischen Systems mit zwei Gleichungen dar.
 */
void euler(double t, double dt, double x[2], dgl RT, dgl BT)
{
  double xtmp[2];
  xtmp[0] = x[0];
  xtmp[1] = x[1];

  x[0] += dt * RT(t, xtmp);
  x[1] += dt * BT(t, xtmp);
}

/**
 * @brief Stellt einen Schritt des Runge-Kutta-Verfahrens vierter Ordnung zur Lösung
 *        des Anfangswertproblems eines dynamischen Systems mit zwei Gleichungen dar.
 */
void runge_kutta4(double t, double dt, double x[2], dgl RT, dgl BT)
{
  double xtmp[2];

  xtmp[0] = x[0];
  xtmp[1] = x[1];
  double f1Rt = RT(t, xtmp);
  double f1Bt = BT(t, xtmp);

  xtmp[0] += dt * 0.5 * f1Rt;
  xtmp[1] += dt * 0.5 * f1Rt;
  double f2Rt = RT(t + dt * 0.5, xtmp);
  double f2Bt = BT(t + dt * 0.5, xtmp);

  xtmp[0] += dt * 0.5 * f2Rt;
  xtmp[1] += dt * 0.5 * f2Rt;
  double f3Rt = RT(t + dt * 0.5, xtmp);
  double f3Bt = BT(t + dt * 0.5, xtmp);

  xtmp[0] += dt * f3Rt;
  xtmp[1] += dt * f3Rt;
  double f4Rt = RT(t + dt, xtmp);
  double f4Bt = BT(t + dt, xtmp);

  x[0] += (dt / 6) * (f1Rt + 2 * f2Rt + 2 * f3Rt + f4Rt);
  x[1] += (dt / 6) * (f1Bt + 2 * f2Bt + 2 * f3Bt + f4Bt);
}

int main(int argc, char *argv[])
{
  gp = openGnuPlot();

  double t_ende = 1000;
  double x_start[2];
  x_start[0] = 30;
  x_start[1] = 80;

  /* solve(0.1, t_ende, x_start, raeuber, beute, euler, "euler_dt01.dat");
  solve(0.01, t_ende, x_start, raeuber, beute, euler, "euler_dt001.dat");
  solve(0.001, t_ende, x_start, raeuber, beute, euler, "euler_dt0001.dat"); */

  /* solve(0.1,   t_ende, x_start, raeuber, beute, runge_kutta2, "RK2_dt01.dat");
  solve(0.01,  t_ende, x_start, raeuber, beute, runge_kutta2, "RK2_dt001.dat");
  solve(0.001, t_ende, x_start, raeuber, beute, runge_kutta2, "RK2_dt0001.dat"); */

  solve(0.1, t_ende, x_start, raeuber, beute, runge_kutta4, "RK4_dt01.dat");
  solve(0.01, t_ende, x_start, raeuber, beute, runge_kutta4, "RK4_dt001.dat");
  solve(0.001, t_ende, x_start, raeuber, beute, runge_kutta4, "RK4_dt0001.dat");

  closeGnuPlot(gp);

  return 0;
}

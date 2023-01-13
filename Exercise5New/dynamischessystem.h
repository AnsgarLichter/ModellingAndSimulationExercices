#include <stdio.h>
#include <stdlib.h>
#include <math.h>


FILE* gp;
/**
 * @brief Öffnet GnuPlot und gibt das Handle zur Pipe zurück.
 */

FILE* openGnuPlot() {
  FILE *gp = popen("gnuplot -p", "w");

  fprintf(gp, "set xlabel 'Räuber';\n");
  fprintf(gp, "set ylabel 'Beute';\n");

  return gp;
}

/**
 * @brief Öffnet einen Plot in GnuPlot.
 *
 * @param plot Das Handle der GnuPlot Pipe
 * @param plotfile Der Dateipfad des Plots
 */

void addPlot(FILE* gp, const char* plotfile) {
  static int isfirst = 1;
  if (isfirst) {
    isfirst = 0;
    fprintf(gp, "plot ");
  } else {
    fprintf(gp, ", ");
  }
  fprintf(gp, "'%s' with linespoints", plotfile);
//   fprintf(gp, "'%s' ", plotfile);
}

/**
 * @brief Schließt die Pipe zu GnuPlot.
 *
 * @param plot Das Handle der GnuPlot Pipe
 */

void closeGnuPlot(FILE* gp) {
  fprintf(gp, ";\n");
  fflush(gp);
  pclose(gp);
}

/**
 * @brief Plottet eine Koordinate in eine GnuPlot Datendatei.
 *
 * @param file Das Handle der GnuPlot Datendatei
 * @param x Die zu plottende X-Koordinate
 * @param y Die zu plottende Y-Koordinate
 */

void plotCoordinate(FILE* file, double x, double y, double t) {
  fprintf(file, "%lf %lf %lf\n", x, y, t);
}

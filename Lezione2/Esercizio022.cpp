#include "error.h"
#include "random.h"
#include "randomwalk.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main(int argc, char *argv[]) {

  int M = 10000; // number of throws
  int n = 100;   // number of blocks
  int L = M / n; // length of the blocks

  int dim = 3;        // number of dimension
  int lengthstep = 1; // step length
  randomwalk RW(dim, lengthstep);
  int nstepmax = 100; // maximum number of step

  // Discrete case
  vector<double> disc_meanRWdistance(n, 0.0); // mean of each block
  double disc_blockmeanRWdistance = 0.0;      // mean of all blocks
  double disc_RWdistance_err = 0.0;           // error of all blocks
  double mean = 0.0;                          // support variable for calculating mean

  ofstream outfilemeanRW("discretemeanRW.dat");
  for (int i = 1; i <= nstepmax; i++) { // loop over number of steps
    disc_blockmeanRWdistance = 0.0;     // loop over blocks
    disc_RWdistance_err = 0.0;
    for (int j = 0; j < n; j++) {
      mean = 0.0;
      for (int k = 0; k < L; k++) {
        RW.restart();                     // reset position in the origin
        RW.makediscretestep(i);           // make i discrete steps
        mean += RW.distancesquared() / L; // calculate mean of one block
      }
      disc_meanRWdistance[j] = sqrt(mean);
      // calculate mean of all blocks
      disc_blockmeanRWdistance += disc_meanRWdistance[j] / n;
    }
    disc_RWdistance_err =  calculate_error(disc_meanRWdistance, n - 1); // calculate error
    // print
    outfilemeanRW << i << "  " << disc_blockmeanRWdistance << "  "  << disc_RWdistance_err << endl;
  }
  outfilemeanRW.close();

  // Continuous case
  // Analogously as above
  vector<double> cont_meanRWdistance(n, 0.0);
  double cont_blockmeanRWdistance = 0.0;
  double cont_RWdistance_err = 0.0;

  ofstream outfilecontmeanRW("continuousmeanRW.dat");
  for (int i = 1; i <= nstepmax; i++) {
    cont_blockmeanRWdistance = 0.0;
    cont_RWdistance_err = 0.0;
    for (int j = 0; j < n; j++) {
      mean = 0.0;
      for (int k = 0; k < L; k++) {
        RW.restart();
        RW.makecontinuousstep(i);
        mean += RW.distancesquared() / L;
      }
      cont_meanRWdistance[j] = sqrt(mean);
      cont_blockmeanRWdistance += cont_meanRWdistance[j] / n;
    }
    cont_RWdistance_err = calculate_error(cont_meanRWdistance, n - 1);
    outfilecontmeanRW << i << "  " << cont_blockmeanRWdistance << "  "  << cont_RWdistance_err << endl;
  }
  outfilecontmeanRW.close();

  RW.RWsaveseed();
  return 0;
}

#include "error.h"
#include "random.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main(int argc, char *argv[]) {
  double S0 = 100;     // initial asset price
  double T = 1;        // delivery time
  double K = 100;      // strike price
  double r = 0.1;      // risk-free interest rate
  double sigma = 0.25; // volatility

  int M = 100000; // number of throws
  int n = 100;    // number of blocks
  int L = M / n;  // block length

  double x = 0.0;                     // support random variable
  double St = 0.0;                    // asset price at t=T.
  vector<double> meanCt(n, 0.0);      // call price mean of each block
  vector<double> meanPt(n, 0.0);      // put price mean of each block
  vector<double> blockmeanCt(n, 0.0); // call price mean of the first i+1 blocks
  vector<double> blockmeanPt(n, 0.0); // put price mean of the first i+1 blocks
  vector<double> Ct_err(n, 0.0);      // call price mean error of the first i+1 blocks
  vector<double> Pt_err(n, 0.0);      // put price mean error of the first i+1 blocks

  // Discrete case
  int t_steps = 100;                      // time steps
  vector<double> St_dt(t_steps + 1, 0.0); // asset price at each time step

  // same vectors as above, but for discrete case
  vector<double> discretemeanCt(n, 0.0);
  vector<double> discretemeanPt(n, 0.0);
  vector<double> discreteblockmeanCt(n, 0.0);
  vector<double> discreteblockmeanPt(n, 0.0);
  vector<double> discreteCt_err(n, 0.0);
  vector<double> discretePt_err(n, 0.0);

  Random rnd;
  int seed[4];
  int p1, p2;
  ifstream Primes("Primes");
  if (Primes.is_open()) {
    Primes >> p1 >> p2;
  } else
    cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();

  ifstream input("seed.in");
  string property;
  if (input.is_open()) {
    while (!input.eof()) {
      input >> property;
      if (property == "RANDOMSEED") {
        input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        rnd.SetRandom(seed, p1, p2);
      }
    }
    input.close();
  } else
    cerr << "PROBLEM: Unable to open seed.in" << endl;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < L; j++) {
      St_dt[0] = S0;
      for (int k = 1; k <= t_steps; k++) {
        x = rnd.Gauss(0, 1); // Generate Gaussian random variable 
        //  Calculate St_dt[k] using the European formula
        St_dt[k] = St_dt[k - 1] * exp((r - sigma * sigma / 2) * T / t_steps +  sigma * x * sqrt(T / t_steps));
      }
      // Calculate the discrete mean call option price and put option price and  add them to the corresponding arrays
      discretemeanCt[i] += exp(-r * T) * (St_dt[t_steps] - K > 0 ? St_dt[t_steps] - K : 0) / L;
      discretemeanPt[i] += exp(-r * T) * (K - St_dt[t_steps] > 0 ? K - St_dt[t_steps] : 0) / L;

      // Calculate St(T) using the European formula
      x = rnd.Gauss(0, 1);
      St = S0 * exp((r - sigma * sigma / 2) * T + sigma * x * sqrt(T));
      // Calculate the mean call option price and put option price and add them to the corresponding arrays
      meanCt[i] += exp(-r * T) * (St - K > 0 ? St - K : 0) / L;
      meanPt[i] += exp(-r * T) * (K - St > 0 ? K - St : 0) / L;
    }

    // calculate the mean and variance of first i+1 blocks
    if (i == 0) {
      blockmeanCt[i] = meanCt[i];
      blockmeanPt[i] = meanPt[i];

      discreteblockmeanCt[i] = discretemeanCt[i];
      discreteblockmeanPt[i] = discretemeanPt[i];
    } else {
      blockmeanCt[i] = blockmeanCt[i - 1] * i / (i + 1) + meanCt[i] / (i + 1);
      blockmeanPt[i] = blockmeanPt[i - 1] * i / (i + 1) + meanPt[i] / (i + 1);

      discreteblockmeanCt[i] = discreteblockmeanCt[i - 1] * i / (i + 1) + discretemeanCt[i] / (i + 1);
      discreteblockmeanPt[i] = discreteblockmeanPt[i - 1] * i / (i + 1) + discretemeanPt[i] / (i + 1);
    }

    // calculate errors
    Ct_err[i] = calculate_error(meanCt, i);
    Pt_err[i] = calculate_error(meanPt, i);

    discreteCt_err[i] = calculate_error(discretemeanCt, i);
    discretePt_err[i] = calculate_error(discretemeanPt, i);
  }

  // print
  ofstream outfilemeanCt("meanCt.dat");
  ofstream outfilemeanPt("meanPt.dat");
  ofstream outfilediscretemeanCt("discretemeanCt.dat");
  ofstream outfilediscretemeanPt("discretemeanPt.dat");

  for (int i = 0; i < n; i++) {
    outfilemeanCt << blockmeanCt[i] << "  " << Ct_err[i] << endl;
    outfilemeanPt << blockmeanPt[i] << "  " << Pt_err[i] << endl;

    outfilediscretemeanCt << discreteblockmeanCt[i] << "  " << discreteCt_err[i]  << endl;
    outfilediscretemeanPt << discreteblockmeanPt[i] << "  " << discretePt_err[i]  << endl;
  }
  outfilemeanCt.close();
  outfilemeanPt.close();
  outfilediscretemeanCt.close();
  outfilediscretemeanPt.close();

  rnd.SaveSeed();
  return 0;
}
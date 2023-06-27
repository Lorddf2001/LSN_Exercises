#include "error.h"
#include "random.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main(int argc, char *argv[]) {

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

  // Exercise 01.1

  int M = 10000; // Number of throws
  int n = 100;   // Number of blocks
  int L = M / n; // Number of throws in each block
  double x;      // Random variable

  vector<double> mean(n, 0.0);          // mean of each block
  vector<double> variance(n, 0.0);      // variance of each block
  vector<double> blockmean(n, 0.0);     // mean of the first i+1 blocks
  vector<double> blockvariance(n, 0.0); // variance of the first i+1 block
  vector<double> blockmeanErr(n, 0.0); // error of the mean of the first i+1 block
  vector<double> blockvarianceErr(n, 0.0); // error of the variance of the first i+1 blocks

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < L; j++) {
      x = rnd.Rannyu();
      mean[i] = mean[i] + x / L;                             // mean block i
      variance[i] = variance[i] + (x - 0.5) * (x - 0.5) / L; // variance block i
    }
    // calculate the mean and variance of first i+1 blocks
    if (i == 0) {
      blockmean[i] = mean[i];
      blockvariance[i] = variance[i];
    } else {
      blockmean[i] = blockmean[i - 1] * i / (i + 1) + mean[i] / (i + 1);
      blockvariance[i] = blockvariance[i - 1] * i / (i + 1) + variance[i] / (i + 1);
    }
    // calculate errors
    blockmeanErr[i] = calculate_error(mean, i);
    blockvarianceErr[i] = calculate_error(variance, i);
  }

  // print
  ofstream outfilemean("mean.dat");
  ofstream outfilevariance("variance.dat");
  for (int i = 0; i < blockmean.size(); i++) {
    outfilemean << blockmean[i] << "  " << blockmeanErr[i] << endl;
    outfilevariance << blockvariance[i] << "  " << blockvarianceErr[i] << endl;
  }
  outfilemean.close();
  outfilevariance.close();

  // Chi squared test
  int M2 = 10000;                 // Random generated numbers
  int n2 = 100;                   // Bins of the histogramm
  vector<double> counts(n2, 0.0); // counts in each bin
  vector<double> Chi2(n2, 0.0);   // Chi2
  int k = 0;                      // random integer value
  int expected = M2 / n2;         // expected value

  for (int i = 0; i < n2; i++) {
    fill(counts.begin(), counts.end(), 0); // reset counts

    for (int j = 0; j < M2; j++) {
      k = rnd.Rannyu() * n2; // random bin index
      counts[k]++;            // increment count
    }
    // calculate chi squared for histogramm
    for (int j = 0; j < n2; j++) {
      Chi2[i] = Chi2[i] + (counts[j] - expected) * (counts[j] - expected) / expected;
    }
  }

  // print
  ofstream outfilechi2("Chi2.dat");
  for (int i = 0; i < Chi2.size(); i++) {
    outfilechi2 << Chi2[i] << endl;
  }

  outfilechi2.close();

  rnd.SaveSeed();
  return 0;
}
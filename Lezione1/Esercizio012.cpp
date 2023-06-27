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

  int M; // # throws
  M = 10000;
  double x;
  int Sum[4] = {1, 2, 10, 100};
  double sum = 0.0;

  // Standard
  ofstream outfilestandard("standard.dat");
  for (int j = 0; j < 4; j++) {
    for (int i = 0; i < M; i++) {
      sum = 0.0;
      for (int k = 0; k < Sum[j]; k++) {
        x = rnd.Rannyu();
        sum = sum + x;
      }
      outfilestandard << sum / Sum[j] << endl;
    }
    outfilestandard << endl;
  }

  outfilestandard.close();

  // Exponential
  double lambda;
  lambda = 1.0; // exponential lambda
  ofstream outfileexp("exp.dat");
  for (int j = 0; j < 4; j++) {
    for (int i = 0; i < M; i++) {
      sum = 0.0;
      for (int k = 0; k < Sum[j]; k++) {
        x = rnd.Exp(lambda);
        sum = sum + x;
      }
      outfileexp << sum / Sum[j] << endl;
    }
    outfileexp << endl;
  }

  outfileexp.close();

  // Lorentz

  double gamma;
  gamma = 1.0;
  double mu;
  mu = 0.0;
  ofstream outfilelorentz("lorentz.dat");
  for (int j = 0; j < 4; j++) {
    for (int i = 0; i < M; i++) {
      sum = 0.0;
      for (int k = 0; k < Sum[j]; k++) {
        x = rnd.Lorentz(gamma, mu);
        sum = sum + x;
      }
      outfilelorentz << sum / Sum[j] << endl;
    }
    outfilelorentz << endl;
  }

  outfilelorentz.close();

  rnd.SaveSeed();
  return 0;
}
#include "error.h"
#include "random.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

// Generator function of Cos(Theta) with Theta in [0, M_PI/2], without using M_PI
double Cos(Random rnd) {
  double x = rnd.Rannyu();        // random x value in [0, 1)
  double y = rnd.Rannyu();        // random y value in [0, 1)
  double r = sqrt(x * x + y * y); // distance from the origin (0,0)
  
  // Reject if the point is outside the unit circle, call the "Cos" function recursively until a valid point is generated.
  // If r=0, the division x/r would result in an error, reject and call the "Cos" function recursively;
  if (r >= 1 || r == 0) {
    return Cos(rnd);
  }
  
  // Uniform random number couples fill uniformly the square, thus the points fallen inside the quarter of the circle (r=1) fill it uniformly and the
  // lines which connects the origin with these points individuate a random angle in [0,M_PI/2]. To evaluate cos(theta) just calculate x/r;
  return x / r;
}

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

  double d = 1;                   // spacing
  double l = 0.5;                 // needle lenght
  int M = 1000000;                // number of throws
  int n = 100;                    // number of blocks
  int L = M / n;                  // number of throws in each block
  int Nhit = 0;                   // number of hits
  vector<double> Pi(n, 0.0);      // approximation of pi in each block
  vector<double> blockPi(n, 0.0); // approximation of pi of the first i+1 blocks
  vector<double> blockPiErr(n, 0.0); // standard deviation of calculated pi
  double endpoint;                   // needle's endpoint
  double height;                     // needle's projection on the y axis

  for (int i = 0; i < n; i++) {
    Nhit = 0; // reset number of hits
    for (int j = 0; j < L; j++) {
      endpoint = rnd.Rannyu();
      height = abs(l * Cos(rnd));
      if (endpoint + height >= d) { // if the needle hits the line
        Nhit++;                     // increase number of hits
      }
    }

    Pi[i] = 2 * l * L / (Nhit * d); // approximation of pi using the probability of hitting a line. P = Nhit/Nthrows, pi = (2 * l) / (P * d)
    // calculate approximation of pi for the first i+1 blocks
    if (i == 0) {
      blockPi[0] = Pi[0];
    } else {
      blockPi[i] = (blockPi[i - 1] * i + Pi[i]) / (i + 1);
    }
    // calculate standard deviation
    blockPiErr[i] = calculate_error(Pi, i);
  }

  // print
  ofstream outputfilePi("Pi.dat");
  for (int i = 0; i < n; i++) {
    outputfilePi << blockPi[i] << "  " << blockPiErr[i] << endl;
  }

  outputfilePi.close();

  rnd.SaveSeed();
  return 0;
}
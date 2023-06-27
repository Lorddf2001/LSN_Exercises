#ifndef __integralMC__
#define __integralMC__

#include "distribution.h"
#include "function.h"
#include "random.h"
#include <cmath>
#include <fstream>
#include <iostream>
using namespace std;

class integralMC {

public:
  integralMC() {
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
    } else {
      cerr << "PROBLEM: Unable to open seed.in" << endl;
    }
  }
  
  ~integralMC() { ; }

  void IntegraSaveSeed(){
    rnd.SaveSeed();
  }

  double IntegraAVE(double xmin, double xmax, function &f, int Ndots) {
    double sum = 0;
    for (int i = 0; i < Ndots; i++) {
      sum = sum + f.Eval(rnd.Rannyu(xmin, xmax));
    }
    return (xmax - xmin) * (sum / Ndots);
  }

  // This function calculates the integral of a function f(x) with respect to a given probability distribution d(x). 
  // The integral is approximated by calculating the average of f(x)/d(x) over Ndots randomly generated points from d(x).
  double IntegraAVEdistribution(function &f, distribution &d, int Ndots) {
    double sum = 0;
    double x;
    for (int i = 0; i < Ndots; i++) {
      x = d.Generator(); // Generate a random point x from the distribution d(x)
      sum = sum + f.Eval(x) / d.Eval(x);
    }
    return (sum / Ndots);
  }

  double IntegraHoM(double xmin, double xmax, double fmax, function &f,  int Ndots) {
    int NHit = 0;
    for (int i = 0; i < Ndots; i++) {
      if (rnd.Rannyu(0, fmax) < f.Eval(rnd.Rannyu(xmin, xmax))) {
        NHit = NHit + 1;
      }
    }
    return (xmax - xmin) * fmax * (double)NHit / Ndots;
  }

private:
  Random rnd;
};

#endif
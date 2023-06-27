#ifndef __distribution__
#define __distribution__

#include "random.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

class distribution {

public:
  virtual double
  Eval(double) const = 0; // pure virtual method for evaluating the distribution
  virtual double Generator() = 0; // pure virtual method for generating a random number form given distribution

  distribution() {
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
  virtual ~distribution() { ; };

protected:
  Random rnd; // random number generator object
};

class linear : public distribution {

public:
  linear() : distribution() { ; }
  ~linear() {}
  double Eval(double x) const override {
    // first order Taylor series expanction of function "fun", normalized;
    return (8 / (M_PI * M_PI)) * (M_PI * M_PI / 4 - M_PI * M_PI / 4 * x);
  }
  double Generator() override {
    // Use the inverse cumulative function to map the random number [0,1) to a point in the distribution (also in [0,1) ).
    return 1 - sqrt(1 - rnd.Rannyu());
  }
};

#endif // __distribution__
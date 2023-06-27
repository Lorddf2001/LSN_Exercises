#include "distribution.h"
#include "error.h"
#include "function.h"
#include "integralMC.h"
#include "random.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int main() {

  fun f;                   // function
  linear lin;              // linear distribution
  integralMC calcintegral; // Monte Carlo integrator

  int M = 10000;                                // number of throws
  int n = 100;                                  // number of blocks
  vector<double> Integral(n, 0.0);              // integral of each block
  vector<double> blockIntegral(n, 0.0);         // integral mean of the first i+1 blocks
  vector<double> IntegralErr(n, 0.0);           // error of the integral mean of the first i+1 blocks
  vector<double> Integralsampling(n, 0.0);      // integral of each block using linear sampling method
  vector<double> blockIntegralsampling(n, 0.0); // integral mean of the first i+1 blocks using linear sampling method
  vector<double> IntegralsamplingErr(n, 0.0);   // error of the integral mean of the first i+1 blocks using linear sampling method

  for (int i = 0; i < n; i++) {
    Integral[i] = calcintegral.IntegraAVE(0, 1, f, M); // calculate integral with M throws

    // calculate the integral mean of first i+1 blocks
    if (i == 0) {
      blockIntegral[i] = Integral[i];
    } else {
      blockIntegral[i] =  blockIntegral[i - 1] * i / (i + 1) + Integral[i] / (i + 1);
    }
    // calculate errors
    IntegralErr[i] = calculate_error(Integral, i);
  }

  // print
  ofstream outfileintegral("Integral.dat");
  for (int i = 0; i < n; i++) {
    outfileintegral << blockIntegral[i] << "  " << IntegralErr[i] << endl;
  }
  outfileintegral.close();

  // sampling linear method
  for (int i = 0; i < n; i++) {
    Integralsampling[i] = calcintegral.IntegraAVEdistribution(f, lin, M); // calculate integral with M throws

    // calculate the integral mean of first i+1 blocks
    if (i == 0) {
      blockIntegralsampling[i] = Integralsampling[i];
    } else {
      blockIntegralsampling[i] = blockIntegralsampling[i - 1] * i / (i + 1) + Integralsampling[i] / (i + 1);
    }
    // calculate errors
    IntegralsamplingErr[i] = calculate_error(Integralsampling, i);
  }

  // print
  ofstream outfileintegralsampling("Integralsampling.dat");
  for (int i = 0; i < n; i++) {
    outfileintegralsampling << blockIntegralsampling[i] << "  "  << IntegralsamplingErr[i] << endl;
  }
  outfileintegralsampling.close();

  calcintegral.IntegraSaveSeed();

  return 0;
}
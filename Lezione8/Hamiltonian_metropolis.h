#ifndef __Hamiltonian_metropolis__
#define __Hamiltonian_metropolis__

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "random.h"


using namespace std;

class Hamiltonian_metropolis{

public: 

  Hamiltonian_metropolis(double mu, double sigma, double delta, int nblk, int nstep){
     ifstream Primes, Seed;

    //Read seed for random numbers
    int p1, p2;
    Primes.open("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();

    Seed.open("seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    Seed.close();

    m_mu = mu;
    m_sigma = sigma;
    m_delta = delta;
    m_nblk = nblk;
    m_nstep = nstep;
    
    x = rnd.Gauss(0, 1);

    return;
  }

  ~Hamiltonian_metropolis() {}

  vector<double> Eigenvalue(int Print)
  {
    vector<double>  H_and_err;    // Vector to store the computed eigenvalue and its error

    for(int iblk=1; iblk <= m_nblk; iblk++) 
    {
      Reset(iblk);   //Reset block averages
      for(int istep=1; istep <= m_nstep; istep++)
      {
        Move();
        Measure();
        Accumulate(); //Update block averages
      }
    H_and_err = Averages(iblk, Print);   //results for current block
    }

    return H_and_err;   // Return the computed eigenvalue and its error
  }

  void Psi2_histo(int nbin, double xmax, int nstep)
  {
    vector<double>  Psi2(nbin+1);     //vector to store bin counts
    double bin_size = 2.0*xmax/nbin;  //Interval [-xmax, xmax]
    int out = 0;                      //support variable to see if the interval [-xmax, xmax] is too little, counts points outside this interval

    for(int i =1; i <= nstep; i ++){
      Move();
      if(x > - xmax || x < xmax){
        Psi2[static_cast<int>((x + xmax) / bin_size)] += 1;   //fill histo
      } else {
        out ++;
      }
    }
    cout << "Number of points out of [ " << - xmax <<", " << xmax << " ] in function Psi2_histo = " << out << endl;

  //Normalization 
    double sum = 0.0;
    for(int i = 0; i <= nbin; i++){
      sum += Psi2[i];
    }
    double mean = sum/nbin;
    for(int i = 0; i <= nbin; i++){
      Psi2[i] = Psi2[i]/(2.0 * xmax * mean);
    }

  //Print 
    ofstream Psi2_hist;
    Psi2_hist.open("Psi2_Histogram.dat");
    for(int i=0; i<=nbin; i++){
      Psi2_hist << -xmax + bin_size * i << "  " << Psi2[i] << endl;
    }
    Psi2_hist.close();
    }

  void SetParameters(double mu, double sigma){
    m_mu = mu;
    m_sigma = sigma;
  }

  void SetBlockParameters(int nblk, int nstep){
    m_nblk = nblk;
    m_nstep = nstep;
  }

  double Unif(){
    return rnd.Rannyu();
  }


  void Move()
  {
      double p;
      double xnew;

      xnew = x + m_delta*(rnd.Rannyu() - 0.5);
      p = Psi2(xnew)/Psi2(x);     //acceptance probability
      if(rnd.Rannyu() < p){
        x = xnew;
        accepted ++;
      }
      attempted ++;
      return;
  }

  double Psi2(double x) //Squared wave function Psi^2
  {
    return pow( exp( - (x-m_mu)*(x-m_mu)/(2.0 * m_sigma*m_sigma) ) + exp( - (x+m_mu)*(x+m_mu)/(2.0 * m_sigma*m_sigma) ), 2);
  }

  void Measure() //Properties measurement
  {
    double Epot = pow(x, 4) - 2.5 * pow(x, 2);
    double Ekin = exp( -(x-m_mu)*(x-m_mu)/(2.0* m_sigma*m_sigma) ) * ( (x-m_mu)*(x-m_mu)/(m_sigma * m_sigma) - 1 )/(m_sigma*m_sigma)
                  + exp( -(x+m_mu)*(x+m_mu)/(2.0* m_sigma*m_sigma) ) * ( (x+m_mu)*(x+m_mu)/(m_sigma * m_sigma) - 1 )/(m_sigma*m_sigma);
    Ekin = - 0.5 * hbar*hbar/m * Ekin /( exp( - (x-m_mu)*(x-m_mu)/(2.0 * m_sigma*m_sigma) ) + exp( - (x+m_mu)*(x+m_mu)/(2.0 * m_sigma*m_sigma) ) );    //hbar, m = 1
   
    // Total energy calculation
    energy = (Epot + Ekin);

    return;
  }

  void Reset(int iblk) //Reset block averages
  {
    if(iblk == 1)
        {
        glob_av = 0;
        glob_av2 = 0;
        }

    blk_av = 0;
   
    blk_norm = 0;
    attempted = 0;
    accepted = 0;
  }


  void Accumulate(void) //Update block averages
  {
    blk_av = blk_av + energy;  
    blk_norm = blk_norm + 1.0;
  }

  vector<double> Averages(int iblk, int Print) //Print results for current block
  {  
    estimate_energy = blk_av/blk_norm; //Energy eigenvalue
    glob_av += estimate_energy;
    glob_av2 += estimate_energy*estimate_energy;
    err_energy = Error(glob_av, glob_av2, iblk);

    if(Print){
      ofstream Metro;
      Metro.open("output_metro.dat",ios::app);
      Metro << iblk << "  " << estimate_energy << "  " << glob_av/(double)iblk << "   "  << err_energy << endl;
      Metro.close();
    }
    
    return {glob_av/(double)iblk, err_energy};
  }

  double Error(double sum, double sum2, int iblk)
  {
    return sqrt(fabs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
  }

// Class members
private: 

// averages
double blk_av, blk_norm, accepted, attempted;
double glob_av, glob_av2;
double estimate_energy;
double err_energy;

// simulation
int m_nstep, m_nblk;
double m_delta;

//configuration
double x;
double xnew;

//parameters, observables
double m_mu, m_sigma;
double hbar = 1;
double m = 1;
double energy;

//Random numbers
int seed[4];
Random rnd;

};

#endif // __Hamiltonian_metropolis__
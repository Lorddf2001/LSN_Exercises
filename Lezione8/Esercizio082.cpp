#include "Hamiltonian_metropolis.h"

int main(){
    double mu_old, sigma_old, mu, sigma;
    int nblk_MC, nstep_MC;
    double delta_MC, alfa, T, Tmin;
    double p;
    vector <double> Eigenvalue_old, Eigenvalue_new, Eigenvalue;
    int Tstep = 0;
    int accepted = 0;
    int attempted = 0;
    ofstream Energy_SA;
    ofstream Parameters_SA;

//----------------------------------------------------------------------------------------------------------------------
// Input
    ifstream ReadInput;
    ReadInput.open("input.in");
    ReadInput >> mu_old;
    ReadInput >> sigma_old;
    ReadInput >> delta_MC;     //energy perturbation for Move() method 
    ReadInput >> nblk_MC;
    ReadInput >> nstep_MC;

    ReadInput >> alfa;     //cooling rate
    ReadInput >> T;
    ReadInput >> Tmin;

    cout << "SA algorithm" << endl;
    cout << "Delta_MC = " << delta_MC << endl;
    cout << "Number of blocks in internal metropolis = " << nblk_MC << endl;
    cout << "Number of steps in each block = " << nstep_MC << endl;

    cout << "Cooling rate Alfa = " << alfa << endl;
    cout << "Initial temperature = " << T << endl;
    cout << "Minimum temperature = " << Tmin << endl;
    
//-------------------------------------------------------------------------------------------------------------------------

    
    Hamiltonian_metropolis H(mu_old, sigma_old, delta_MC, nblk_MC, nstep_MC);
    Energy_SA.open("output_Energy_SA.dat");
    Parameters_SA.open("output_Parameters_SA.dat");

    while(T > Tmin){
        Eigenvalue_old = H.Eigenvalue(0);               // Compute the eigenvalue

        mu = mu_old + 0.1 * (H.Unif() - 0.5);           // Update mu using a random perturbation
        sigma = sigma_old + 0.1 * (H.Unif() - 0.5);     // Update sigma using a random perturbation
        H.SetParameters(mu, sigma);                     // Set the new values of mu and sigma 

        Eigenvalue_new = H.Eigenvalue(0);               // Compute the eigenvalue with the updated parameters
        p = exp( - 1./T * (Eigenvalue_new[0] - Eigenvalue_old[0]));     // Calculate the acceptance probability

        if(H.Unif() < p){   // Accept the new parameters with probability p
            Eigenvalue = Eigenvalue_new;
            mu_old = mu;
            sigma_old = sigma;
            accepted ++;
        } else {    // Reject the new parameters and revert to the old values
            mu = mu_old;
            sigma = sigma_old;
            H.SetParameters(mu, sigma);
            Eigenvalue = Eigenvalue_old;
        }
        attempted ++;
        Tstep ++;
        T = alfa * T;       // Decrease the temperature according to the cooling schedule
        //Print
        Energy_SA << Eigenvalue[0] << "   " << Eigenvalue[1] << endl;
        Parameters_SA << mu << "   " << sigma << endl;
    }

    Energy_SA.close();
    Parameters_SA.close();
    
    cout << "Temperature steps " << Tstep << endl;
    cout << endl;
    cout << "Acceptance rate  " << (double)accepted/attempted << endl;
    cout << "Best values  " << mu << "  " << sigma << endl;
    cout << "Best energy  " << Eigenvalue[0] << "  " <<  Eigenvalue[1] << endl;

    Eigenvalue = H.Eigenvalue(1);   // Compute and print the final eigenvalue

    H.Psi2_histo(120, 3, 100000);   // Generate a histogram of the wave function
    return 0;
}
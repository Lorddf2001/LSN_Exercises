#include "map.h"
#include "route.h"
#include "population.h"

int main(){
    //Random 
    Random rnd;
    int seed[4];
    ifstream Primes, Seed;
    int p1, p2;
    Primes.open("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();

    Seed.open("seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    Seed.close();

    //Input 
    ifstream Input;
    Input.open("input.in");
    int n_population, n_city, shape, nstep;
    string Shape;
    Input >> n_population;
    Input >> n_city;
    Input >> shape;
    Input >> nstep;
    Input.close();
    
    if (shape == 0){
        Shape = "circle";
    }else{
        Shape = "square";
    }
    
    ofstream Length;
    Length.open("output_lengthmean_"+Shape+".dat");
    
    //Simulation
    Population pop(n_population, n_city, shape);
    for(int i  = 0; i < nstep; i++){
        pop.nextgen();
        pop.mutation();
        pop.sorting();
        Length << pop.mean_length() << endl;
    }
    Length.close();
    pop.print_map(shape);
    pop.print_bestroute(shape);

    return 0;
}


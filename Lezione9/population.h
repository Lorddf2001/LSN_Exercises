#ifndef __Population__
#define __Population__

#include "map.h"
#include "route.h"

class Population{

public:

    Population(int npopulation, int ncities, int shape){
        //Read seed for random numbers
        ifstream Primes, Seed;
        int p1, p2;
        Primes.open("Primes");
        Primes >> p1 >> p2 ;
        Primes.close();

        Seed.open("seed.in");
        Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        rnd.SetRandom(seed,p1,p2);
        Seed.close();

        //Generate map
        Map temp_map(shape, ncities, rnd);
        m_map = temp_map;

        //Generate routes
        m_npopulation = npopulation;
        m_pop.resize(m_npopulation);
        for(int i = 0; i < m_npopulation; i++){
            Route route_i(ncities, rnd, m_map);
            m_pop[i] = route_i;
        }
    }

    ~Population () {}

    void change_map(int ncities, int shape){
        Map temp_map(shape, ncities, rnd);
        m_map = temp_map;

        for(int i = 0; i < m_npopulation; i++){
            Route route_i(ncities, rnd, m_map);
            m_pop[i] = route_i;
        }
    }

    void nextgen(){
        int father_index, mother_index;
        for(int i = 0; i < m_npopulation/2; i++){        //n/2 couples
            father_index = selection();
//            cout << "selection 1 " << father_index << endl;
            mother_index = selection();
//            cout << "selection 2 " << mother_index << endl;
            crossover(father_index, mother_index);
        }
        sorting();  
        m_pop.erase(m_pop.begin() + m_map.getncity() , m_pop.begin() + 2.0 * m_map.getncity());     // Keep only the best individuals
    }

    int selection(){    //Select an individual from the population based on their fitness level, giving preference to individuals with higher fitness values.
        double p = 4.0;
        return (int)( m_npopulation * pow(rnd.Rannyu(), p));
    }

    void mutation(){
        double mutation;
        for(int i=0; i < m_npopulation; i++ ){
            mutation = rnd.Rannyu();
            if(mutation < 0.01){
                m_pop[i].shift(rnd);
            } else if( mutation < 0.03){
                m_pop[i].inversion(rnd);
            } else if( mutation < 0.04 ){
                m_pop[i].block_swap(rnd);
            } else if( mutation < 0.07){
                m_pop[i].single_swap(rnd);
            }
        }
    }

    void crossover (int father, int mother) {
        // Create children as a copy of the parents
        Route child1 = m_pop[father];
        Route child2 = m_pop[mother];
        // crossover probability
        if( rnd.Rannyu() < 0.90){
            int support1 = 0;
            int support2 = 0;
            int cut = rnd.Discrete(1, m_map.getncity() - 2);    // Randomly select a cutting point

            // Resize children's route up to the cutting point
            vector<int> child1_route;
            vector<int> child2_route;
            child1_route = child1.getroute();
            child1_route.resize(cut);
            child2_route = child2.getroute();
            child2_route.resize(cut);

            // Perform crossover by combining the parent routes up to the cutting point
            for(int i = 0; i < m_map.getncity(); i ++){    // Iterate over cities in the parent routes
                support1 = 0;
                support2 = 0;
                // Check if the city from the mother parent already exists in child1's route
                for(int j = 0; j < child1_route.size() ; j++){
                    if( m_pop[mother].getcityindex(i) == child1_route[j] ){
                        support1 ++;
                        break;
                    }
                }
                // If city doesn't exist, add it to child1's route
                if (support1 == 0 ) child1_route.push_back(m_pop[mother].getcityindex(i));

                // Check if the city from the father parent already exists in child2's route
                for(int j = 0; j < child2_route.size(); j++){
                    if( m_pop[father].getcityindex(i) == child2_route[j] ){
                        support2 ++;
                        break;
                    }
                }
                // If city doesn't exist, add it to child2's route
                if (support2 == 0 ) child2_route.push_back(m_pop[father].getcityindex(i));
            }
            // Update children's route
            child1.updateroute(child1_route);
            child2.updateroute(child2_route);
        }
        // Add the modified child individuals to the population
        m_pop.push_back(child1);
        m_pop.push_back(child2);
    }

    void sorting(){
        for(int i = 0; i < m_pop.size(); i++){
            m_pop[i].recalculate_length(m_map);
        }
        sort(m_pop.begin(), m_pop.end());       //with overridden < 
    }

    double route_length(int iroute){
        return m_pop[iroute].getlength();
    }

    void print_routes(int shape){
        for(int i = 0; i < m_npopulation; i++){
            m_pop[i].print(shape);
        }
    }

    void print_bestroute(int shape){    //require sorting before calling
        m_pop[0].print(shape);
    }

    void print_map(int shape){
        m_map.print_position(shape);
    }

    void print_city_distances(){
        m_map.print_city_distance();
    }

    void print_route_length(){
        for(int i = 0; i < m_npopulation; i++){
            m_pop[i].print_length();
        }
    }

    double mean_length(){   //of half population
        double sum = 0.0;
        for(int i = 0; i < m_npopulation/2; i++){
            sum += m_pop[i].getlength();
        }
        return sum/(m_npopulation/2);
    }

    int getncity(){
        return m_map.getncity();
    }

private:

int m_npopulation;
vector<Route> m_pop;
Map m_map;

int seed[4];
Random rnd;
};

#endif // __Population__
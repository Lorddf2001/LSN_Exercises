#ifndef __Population__
#define __Population__

#include "map.h"
#include "route.h"

class Population{

public:

    Population(int npopulation,Random& rnd){
        m_rnd = rnd;

        Map temp_map;
        m_map = temp_map;

        m_npopulation = npopulation;
        m_pop.resize(m_npopulation);
        for(int i = 0; i < m_npopulation; i++){
            Route route_i(m_rnd, m_map);
            m_pop[i] = route_i;
        }
    }

    ~Population () {}

/*  useless here
    void change_map(int ncities, int shape){
        Map temp_map(shape, ncities, m_rnd);
        m_map = temp_map;

        for(int i = 0; i < m_npopulation; i++){
            Route route_i(ncities, m_rnd, m_map);
            m_pop[i] = route_i;
        }
    }
*/
    void nextgen(){
        int father_index, mother_index;
        for(int i = 0; i < m_npopulation/2; i++){        //n/2 couples
            father_index = selection();
            mother_index = selection();
            crossover(father_index, mother_index);
        }
        sorting();  
        m_pop.erase(m_pop.begin() + m_map.getncity() , m_pop.begin() + 2.0 * m_map.getncity());     //only the best remains
    }

    int selection(){
        double p = 4.0;
        return (int)( m_npopulation * pow(m_rnd.Rannyu(), p));
    }

    void mutation(){
        double mutation;
        for(int i=0; i < m_npopulation; i++ ){
            mutation = m_rnd.Rannyu();
            if(mutation < 0.01){
                m_pop[i].shift(m_rnd);
            } else if( mutation < 0.03){
                m_pop[i].inversion(m_rnd);
            } else if( mutation < 0.04 ){
                m_pop[i].block_swap(m_rnd);
            } else if( mutation < 0.07){
                m_pop[i].single_swap(m_rnd);
            }
        }
    }

    void crossover (int father, int mother) {
        int cut = m_rnd.Discrete(1, m_map.getncity() - 2);

        // Create children as a copy of the parents
        Route child1 = m_pop[father];
        Route child2 = m_pop[mother];
        // Resize child's route to the cutting point
        child1.getroute().resize(cut);
        child2.getroute().resize(cut);
        
        // Perform crossover by combining the parent routes up to the cutting point
        for(int i = 0; i < m_map.getncity(); i ++){    // Add cities from the parent if not already present in child
            if(! child1.check_city(m_pop[father].getcityindex(i)) ) child1.getroute().push_back(m_pop[father].getcityindex(i));     
            if(! child2.check_city(m_pop[mother].getcityindex(i)) ) child2.getroute().push_back(m_pop[mother].getcityindex(i));
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

    void print_routes(){
        for(int i = 0; i < m_npopulation; i++){
            m_pop[i].print();
        }
    }

    void print_bestroute(){
        m_pop[0].print();
    }

    void print_map(){
        m_map.print_position();
    }

    void print_city_distances(){
        m_map.print_city_distance();
    }

    void print_route_length(){
        for(int i = 0; i < m_npopulation; i++){
            m_pop[i].print_length();
        }
    }

    void print_bestroute_lenght(){
        m_pop[0].print_length();
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

//Methods for parallel algorithm

    Route getbest(){
        return m_pop[0];
    }

    void exchangebest(Route best){
        m_pop[m_npopulation - 1] = best;
    }

    void insert(Route route, int i){
        m_pop[i] = route;
    }

private:

int m_npopulation;
vector<Route> m_pop;
Map m_map;

Random m_rnd;
};

#endif // __Population__
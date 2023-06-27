#ifndef __Route__
#define __Route__

#include <algorithm>
#include "map.h"

using namespace std;

class Route{

public:

    Route(){}

    Route(Random& rnd, Map& map){
        m_ncity = 50;
        m_route.resize(m_ncity);  
        
        //Generate random route
        for (int i = 0; i < m_ncity; i++) {
            m_route[i] = i;     //fill route 
        }
        //shuffle
        int j;
        for (int i = m_ncity - 1; i > 0; i--) {     //first city 0  
            j = rnd.Rannyu(1, i);
            swap(m_route[i], m_route[j]);
        }

        //calculate length
        double length = 0.0;
        for(int i = 0; i < m_ncity - 1; i++){
            length += map.getdistance(m_route[i], m_route[i+1]);
        }
        length += map.getdistance(m_route[m_ncity - 1], m_route[0]);
        m_length = length;
    }

    //Constructor for sending_order in parallel algorithm
    Route(Random& rnd, int n){      
        m_ncity = n;
        m_route.resize(m_ncity);
        //Generate random route
        for (int i = 0; i < m_ncity; i++) {
            m_route[i] = i;     //fill route
        }
        //shuffle
        int j;
        for (int i = m_ncity - 1; i > 0; i--) {       
            j = rnd.Rannyu(0, i);           //Attention: the first "city" can be different from 0
            swap(m_route[i], m_route[j]);
        }
        m_length = 0;    //useless for this type of route, fill with 0
    }

    ~Route () {}

    void check(){
        if(0 != m_route[0]){
            cout << "Route doesn't start from city number 0" << endl;
        }
        for(int i = 0; i < m_ncity; i++){
            for(int j = i+1; j < m_ncity; j++){
                if(m_route[i] == m_route[j]){
                    cout << "Same city visited twice" << endl;
                    break;
                }
            }
        }
        
        
    }

    vector<int> getroute(){
        return m_route;
    }

    int getcityindex(int index){
        return m_route[index];
    }

    int getncity(){
        return m_ncity;
    }

    double getlength() const{
        return m_length;
    }

    void recalculate_length(Map& map){
        double length = 0.0;
        for(int i = 0; i < m_ncity - 1; i++){
            length += map.getdistance(m_route[i], m_route[i+1]);
        }
        length += map.getdistance(m_route[m_ncity - 1], m_route[0]);
        m_length = length;
    }

    void single_swap(Random& rnd){
        int i, j;     
        do{
            i = rnd.Discrete(1, m_ncity);
            j = rnd.Discrete(1, m_ncity);
        } while( i == j );
        swap(m_route[i], m_route[j]);
    }

    void shift(Random& rnd){
        int shift_length = rnd.Discrete(1, m_ncity);                  //lenght of the shift
        int city_to_shift = rnd.Discrete(1, m_ncity);                 //number of cities to shift
        int start = rnd.Discrete(1, m_ncity);                         //Position of first city to shift
        for(int i = 0; i < shift_length; i++){
            for(int j = city_to_shift - 1; j >= 0; j--){
                swap(m_route[pbc(start + i + j)], m_route[pbc(start + i + j + 1)]);
            }
        }
    }
           
    void block_swap(Random& rnd){  
      int start =  rnd.Discrete(1, m_ncity);
      int city_to_move = rnd.Discrete(2, m_ncity/2);

      for(int i = 0; i < city_to_move; i++){
        swap(m_route[pbc(start + i)], m_route[pbc(start + i + city_to_move)]);
      }
    }

    void inversion(Random& rnd){
        int city_to_invert = rnd.Discrete(2, m_ncity);
        int start = rnd.Discrete(1, m_ncity);
        for (int i = 0; i < city_to_invert / 2; i++) {
        swap(m_route[pbc(start + i)], m_route[pbc(start + city_to_invert - i - 1)]);
        }
    }

    void print(){
        ofstream Route;
        Route.open("route_.dat", ios::app);
        for(int i=0; i < m_ncity; i++ ){
            Route << m_route[i] << "  ";
        }
        Route << endl;
        Route.close();
    }

    void print_length(){
        ofstream Length;
        Length.open("route_length.dat", ios::app);
        Length << m_length << endl;
        Length.close();
    }

    int pbc(int icity){
        if(icity == 0) return 1;
        while(icity >= m_ncity){
            icity -= (m_ncity - 1);
        }
        return icity;
    }

    bool check_city(int city_index){
        for(int i = 0; i < m_ncity; i++){
            if(m_route[i] == city_index) return true;
        }
        return false;
    }

    bool operator<(const Route& other) const{
        return getlength() < other.getlength();
    }

    //Methods for parallel algorithm
    
    void reshuffle(Random& rnd){
        //shuffle
        int j;
        for (int i = m_ncity - 1; i > 0; i--) {      
        j = rnd.Rannyu(0, i);                   //Attention: the first "city" can be different from 0
        swap(m_route[i], m_route[j]);
        }
    }

    void updateroute(vector<int> route){
        for(int i = 0; i < m_ncity; i++){
            m_route[i] = route[i];
        }
    }

private:
int m_ncity;
vector<int> m_route;
double m_length;

};

#endif // __Route__
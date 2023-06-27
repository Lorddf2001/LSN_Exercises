#ifndef __Map__
#define __Map__

#include <iostream>
#include <fstream>
#include <sstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <string>
#include "random.h"


using namespace std;

class Map {

public: 
    
    Map(){
        ifstream America;
        America.open("American_capitals.in");
        m_ncity = 50; //American capitals
        m_position.resize(2*m_ncity);
        if (America.is_open()) {
            for (int i = 0; i < m_ncity; i++) {
                America >> m_position[2 * i] >> m_position[2 * i + 1];
            }
        America.close();
        } else {
            cout << "The file American_capitals.in cannot be opened." << endl;
        }

        //Calculate distances 
        m_distances.resize(m_ncity, vector<double>(m_ncity));
        double x1, x2, y1, y2, distance;
        for (int i = 0; i < m_ncity; i++) {
            x1 = m_position[2*i];
            y1 = m_position[2*i+1];

            for (int j = i+1; j < m_ncity; j++) {
                x2 = m_position[2*j];
                y2 = m_position[2*j+1];
                distance = sqrt( pow(x2 - x1, 2) + pow(y2 - y1, 2) );

                //Store in m_distances matrix
                m_distances[i][j] = distance;
                m_distances[j][i] = distance; // Copy to the lower half of the matrix
            }
        }
    }

    ~Map () {}

    double getdistance(int i, int j){
        return m_distances[i][j];
    }

    int getncity(){
        return m_ncity;
    }

    vector<double> getposition(){
        return m_position;
    }

    void print_city_distance(){
        ofstream Matrix;
        Matrix.open("distances.dat");
        for(int i = 0; i < m_ncity; i++) {
            for (int j = 0; j < m_ncity; j ++) {
                Matrix << setw(10) << m_distances[i][j] << " ";
            }
        Matrix << endl;
        };
        Matrix.close();
    }

    void print_position(){
        ofstream Position;
        Position.open("map.dat");
            for(int i = 0; i < m_ncity; i++ ){
                Position << m_position[2*i] << "  " << m_position[2*i + 1] << endl;  
            }
        Position.close();
    }

private:
int m_ncity;    //number of cities
vector<double> m_position;  //Position vector
vector<vector<double>> m_distances;     //distances matrix

};


#endif // __Map__
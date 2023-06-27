#ifndef __Map__
#define __Map__

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <string>
#include "random.h"


using namespace std;

class Map {

public: 

    Map(){}
    
    Map(int shape, int ncity, Random& rnd){
        m_ncity = ncity;
        m_position.resize(2 * m_ncity);
        
        if(shape == 0){ //circumference
            double angle;
            for(int i = 0; i < m_ncity; i ++){
                angle = rnd.Rannyu(0, 2.0 *M_PI);
                m_position[2*i] = cos(angle);     //X_i
                m_position[2*i+1] = sin(angle);   //Y_i
            }
        } else {    //square
            for(int i = 0; i < m_ncity; i ++){
                m_position[2*i] = rnd.Rannyu(-1, 1);     //X_i
                m_position[2*i+1] = rnd.Rannyu(-1, 1);   //Y_i
            }
        }

        //Calculate distances squared
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
                m_distances[j][i] = distance; // Copy to the lower half of the matrix, simmetric
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

    void print_position(int shape){
        string Shape;
        if (shape == 0){
        Shape = "circle";
        }else{
        Shape = "square";
        }

        ofstream Position;
        Position.open("map_"+Shape+".dat");
            for(int i = 0; i < m_ncity; i++ ){
                Position << m_position[2*i] << "  " << m_position[2*i + 1] << endl;  
            }
        Position.close();
    }

private:
int m_ncity;                            //Number of cities
vector<double> m_position;              //Position vector
vector<vector<double>> m_distances;     //Distances matrix
};


#endif // __Map__
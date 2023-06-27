#ifndef __error__
#define __error__

#include <vector>
#include <cmath>
#include <iostream>

double calculate_error(const std::vector<double>& values, int N) {
    if ( (int)values.size() <= N) {
        std::cerr << "Error: N is too large." << std::endl;
        return 0;
    }

    if( N == 0 ) { return 0; }

    // Resize the input vector to size N
    std::vector<double> Support(values.begin(), values.begin() + N + 1);

    // Calculate the squared values and the mean
    double mean = 0;
    double squared_mean = 0;
    for(int i = 0; i <= N; i++) {
        mean += Support[i]/(N + 1);
        squared_mean += Support[i] * Support[i]/(N+1);
    }

    // Calculate and return the error
    double error = sqrt((squared_mean - mean * mean) / N);
    return error;
}

#endif // __error__
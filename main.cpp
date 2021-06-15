#include <iostream>
#include <random>
#include <iomanip>
#include <cmath>
#include <map>
#include <ctime>
#include <string>
#include "cubeBased.h"
#include "monteCarlo.h"
std::map<int, double> Volume = {{1, 2.0},
                                {2, 3.14159265358979},
                                {3, 4.18879020478639},
                                {4, 4.93480220054468},
                                {5, 5.26378901391432},
                                {6, 5.16771278004997},
                                {7, 4.7247659703314},
                                {8, 4.05871212641677},
                                {9, 3.29850890273871}};

// double get_actual_volumn(int dimension) {
//    switch (dimension) {
//        case 1: return 2;
//        case 2: return M_PI;
//        case 3: return 4*M_PI/3;
//        case 4: return pow(M_PI,2)/2;
//        case 5: return 8*pow(M_PI,2)/15;
//        case 6: return pow(M_PI,3)/6;
//        case 7: return 16*pow(M_PI,3)/105;
//        case 8: return pow(M_PI,4)/24;
//        case 9: return 32*pow(M_PI,4)/945;
//        default: return 0.0;
//    }
//}


//void generate_random_point(std::mt19937_64, std::uniform_real_distribution<double>,int, double*);
//
//double euclidean_distance(const double* , int);
//
//bool phi(const double*, int);
//
//void monte_carlo_estimate(int , double);
//
//void cube_based_estimate(int, double);

int main(int argc, char** argv) {
    int dimension = 0;
    if (argc > 1) {
        dimension = atoi(argv[1]);
    } else {
        std::cout << "Please specify a dimention as a command line arg." << std::endl;
        return 0;
    }
    if (dimension < 2 or dimension > 9) {
        std::cout << "Dimension " << dimension << "not supported." << std::endl;
        return 0;
    }

    bool fixedN = false;
    std::string s = "fixedN";
    if (argc == 3 && s.compare(argv[2]) == 0) {
        fixedN = true;
    }

    std::cout << "Monte Carlo Estimation: " << std::endl;
    MonteCarloEstimator monteCarlo(dimension, Volume.at(dimension), fixedN);
    monteCarlo.monte_carlo_estimate();

    if (dimension < 3) {
        std::cout << "\nCube Based Estimation: " << std::endl;
        CubeBasedEstimator cubeBased(dimension, Volume.at(dimension), fixedN);
        cubeBased.cube_based_estimate();
    }
    return 0;
}

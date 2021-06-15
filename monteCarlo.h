//
// Created by Yinan Zhou on 6/10/21.
//

#ifndef FINAL_CS206_MONTECARLO_H
#define FINAL_CS206_MONTECARLO_H
#include <random>

class MonteCarloEstimator
{
public:
    int dimension;
    double target_volume;
    std::mt19937_64 generator;
    std::uniform_real_distribution<double> distribution;
    double M;

    MonteCarloEstimator(int,double,bool);
    bool phi();
    void monte_carlo_estimate();
};

#endif //FINAL_CS206_MONTECARLO_H

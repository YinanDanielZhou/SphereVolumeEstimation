//
// Created by Yinan Zhou on 6/9/21.
//

#include <iostream>
#include <random>
#include <iomanip>
#include <cmath>
#include <map>
#include <ctime>
#include "monteCarlo.h"

std::map<int, double> M_size = {{1, 100},
                                {2, 2 * pow(10, 8)},
                                {3, 5 * pow(10, 8)},
                                {4, 1 * pow(10, 9)},
                                {5, 2 * pow(10, 9)},
                                {6, 5 * pow(10, 9)},
                                {7, 1 * pow(10, 10)},
                                {8, 5 * pow(10, 9)},
                                {9, 5 * pow(10, 9)}};

MonteCarloEstimator::MonteCarloEstimator(int d, double v, bool f) {
    dimension = d;
    target_volume = v;
    std::random_device rd;
    std::mt19937_64 g(rd());
    std::uniform_real_distribution<double> dist(-1.0,1.0);

    generator = g;
    distribution = dist;

    M = M_size.at(dimension);

    // the case for fixed N = 10^6
    if (f) {
        double fixed_N = pow(10,6);
        M = fixed_N;
    }

    // for 8 digit precision
//    M = 2 * pow(10, 10);
}

bool MonteCarloEstimator::phi() {
    double sumSquares = 0;
    for (int d = 0; d < dimension; d++) {
        double value = distribution(generator);
        sumSquares += value * value;
        if (sumSquares > 1) {
            return false;
        }
    }
    return true;
}

void MonteCarloEstimator::monte_carlo_estimate() {
    double point[dimension];

    int k = 100;

    double sample_results [k];
    double sample_result_sum = 0;
    double single_sample_result;

    std::clock_t start;
    double duration;
    start = std::clock();

    double sample_mean;
    double sample_variance;
    double final_variance;
    double std_deviation;

    for (int run = 0; run < k; run++) {
        single_sample_result = 0.0;
        for (double n = 0; n < M/k; n++) {
//            generate_random_point(generator, distribution, dimension, point);
            if (phi()) {
                single_sample_result += 1.0;
            }
        }
        single_sample_result = pow(2, dimension) * single_sample_result * k / M;

        sample_results[run] = single_sample_result;
        sample_result_sum += single_sample_result;
//        std::cout << std::setprecision(15) << run << "   "  << single_sample_result << std::endl;

//        sample_mean = sample_result_sum / (run + 1);
////        std::cout << "sample mean: " << std::setprecision(15) << sample_mean << std::endl;
////        std::cout << "accuracy: " << std::setprecision(15) << std::abs(sample_mean - target_volume) << std::endl;
//
//        sample_variance = 0;
//        for (int r = 0; r < run; r++) {
//            sample_variance += pow(sample_results[r] - sample_mean, 2);
//        }
//        if (run > 0) {
//            sample_variance = sample_variance/ run;
//        }
//
//        final_variance = sample_variance / (run + 1);
//        std_deviation = sqrt(final_variance);

//        std::cout << "precision: " << std::setprecision(15) << 2 * std_deviation / sample_mean << std::endl;
    }

    sample_mean = sample_result_sum / k;
    sample_variance = 0;
    for (int r = 0; r < k; r++) {
        sample_variance += pow(sample_results[r] - sample_mean, 2);
    }
    sample_variance = sample_variance/ (k-1);
    final_variance = sample_variance / k;
    std_deviation = sqrt(final_variance);

    std::cout << "sample mean: " << std::setprecision(15) << sample_mean << std::endl;
    std::cout << "accuracy: " << std::setprecision(15) << std::abs(sample_mean - target_volume) << std::endl;
    std::cout << "precision: " << std::setprecision(15) << 2 * std_deviation / sample_mean << std::endl;

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    std::cout<<"execution time: "<< duration <<'\n';
}
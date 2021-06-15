//
// Created by Yinan Zhou on 6/9/21.
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
#include <ctime>
#include <vector>
#include "cubeBased.h"

std::map<int, double> K_size = {{1, 100},
                                {2, 1.5 * pow(2,17)},
                                {3, 1.3 * pow(2,11)},
                                {4, 1.3 * pow(2,11)},
                                {5, 1.3 * pow(2,11)},
                                {6, 1.3 * pow(2,11)},
                                {7, 1.3 * pow(2,11)},
                                {8, 1.3 * pow(2,11)},
                                {9, 1.3 * pow(2,11)},
                                {10, 1.3 * pow(2,11)}};

double log_frequency = pow(10,8);

CubeBasedEstimator::CubeBasedEstimator(int d, double v, bool f) {
    dimension = d;
    target_volume = v;
    K = K_size.at(dimension);

    // the case for fixed N = 10^6
    if (f) {
        double fixed_N = pow(10,6);
        K = std::round(std::pow(fixed_N, 1.0/dimension));
    }

    cube_l = 2.0/K;
    max_to_surface = sqrt(dimension) * cube_l / 2;
    in_threshold = 1 - max_to_surface;
    out_threshold = 1 + max_to_surface;

    for (int d = 0; d < dimension; d++) {
        cube_center.push_back(-1 + cube_l/2);
    }
    original_cube_center = cube_center;

    in_cubes_count = 0;
    out_cubes_count = 0;
    processed_cubes_count = 0;
    total_cubes = pow(K,dimension);
}

bool CubeBasedEstimator::has_vertex_in() {
    // given a cube with center outside the sphere, does it have a vertex in the sphere?
    // find the vertex closest to the origin and check that vertex
    std::vector<double> vertex = cube_center;
    double vertex_distance = 0;
    for (int d = 0; d < dimension; d++) {
        if (vertex.at(d) > 0) {
            vertex.at(d) = vertex.at(d) - cube_l/2;
        } else {
            vertex.at(d) = vertex.at(d) + cube_l/2;
        }
        vertex_distance += vertex.at(d) * vertex.at(d);
        if (vertex_distance > 1) {
            return false;
        }
    }
    return true;
}

bool CubeBasedEstimator::has_vertex_out() {
    // given a cube with center inside the sphere, does it have a vertex outside the sphere?
    // find the vertex furthest from the origin and check that vertex
    std::vector<double> vertex = cube_center;
    double vertex_distance = 0;
    for (int d = 0; d < dimension; d++) {
        if (vertex.at(d) >= 0) {
            vertex.at(d) = vertex.at(d) + cube_l/2;
        } else {
            vertex.at(d) = vertex.at(d) - cube_l/2;
        }
        vertex_distance += vertex.at(d) * vertex.at(d);
        if (vertex_distance > 1) {
            return true;
        }
    }
    return false;
}

int CubeBasedEstimator::phi_cube() {
    double center_distance = 0;
    for (int d = 0; d < dimension; d++) {
        center_distance += cube_center.at(d) * cube_center.at(d);
    }
    center_distance = sqrt(center_distance);

    if (center_distance >= out_threshold) return -1;
    if (center_distance <= in_threshold) return 1;
    if (center_distance > 1 && !has_vertex_in()) return -1;
    if (center_distance < 1 && !has_vertex_out()) return 1;
    return 0;
}

void CubeBasedEstimator::iterate_through_all_cubes(int current_d) {
    int cube_type;
    // reset
    cube_center.at(current_d-1) = original_cube_center.at(current_d-1);
    if (current_d == dimension) {
        int k = 0;
        do {
            cube_type = phi_cube();
            processed_cubes_count += 1;
            if (cube_type == 1) {
                in_cubes_count += 1;
            } else if (cube_type == -1) {
                out_cubes_count += 1;
            }
            cube_center.at(current_d - 1) = cube_center.at(current_d - 1) + cube_l;
            k += 1;
        } while (k < K);
    }
    else {
        for (int k = 0; k < K - 1; k++) {
            iterate_through_all_cubes(current_d + 1);
            cube_center.at(current_d - 1) = cube_center.at(current_d - 1) + cube_l;
        }
        iterate_through_all_cubes(current_d + 1);
    }
//    if (std::fmod(processed_cubes_count, log_frequency) == 0) {
//        std::cout << "processed: " << processed_cubes_count << std::endl;
//    }
}

void CubeBasedEstimator::cube_based_estimate() {
    std::clock_t start;
    double duration;
    start = std::clock();
    std::cout << "total cubes: " << total_cubes << std::endl;

    iterate_through_all_cubes(1);
    double whole_cube_volume = pow(2,dimension);

//    std::cout << in_cubes_count << " " << out_cubes_count << std::endl;

    double on_cubes_count = total_cubes - in_cubes_count - out_cubes_count;
//    std::cout << std::setprecision(15) << "In: " << in_cubes_count/total_cubes<< std::endl;
//    std::cout << std::setprecision(15) << "Out: " << out_cubes_count/total_cubes<< std::endl;
//    std::cout << std::setprecision(15) << "On: " << on_cubes_count/total_cubes<< std::endl;


    double upper_bound = whole_cube_volume * (total_cubes - out_cubes_count) / total_cubes;
    double lower_bound = whole_cube_volume * in_cubes_count / total_cubes;
    double estimate = whole_cube_volume * (in_cubes_count + on_cubes_count/2) / total_cubes;

    std::cout << std::setprecision(15) << "lower bound: " << lower_bound << std::endl;
    std::cout << std::setprecision(15) << "upper bound: " << upper_bound << std::endl;
    std::cout << std::setprecision(15) << "Precision: " << upper_bound - lower_bound << std::endl;
    std::cout << std::setprecision(15) << "Estimate: " << estimate << ", accuracy: " << std::abs(estimate - target_volume) << std::endl;

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    std::cout<<"execution time: "<< duration <<'\n';

}

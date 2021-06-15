//
// Created by Yinan Zhou on 6/10/21.
//

#ifndef FINAL_CS206_CUBEBASED_H
#define FINAL_CS206_CUBEBASED_H

#include <vector>

class CubeBasedEstimator
{
public:
    int dimension;
    double target_volume;
    double K;
    double cube_l;
    double max_to_surface;
    double in_threshold;
    double out_threshold;
    std::vector<double> cube_center;
    std::vector<double> original_cube_center;
    double in_cubes_count;
    double out_cubes_count;
    double processed_cubes_count;
    double total_cubes;

    CubeBasedEstimator(int,double,bool);
    bool has_vertex_in();
    bool has_vertex_out();
    int phi_cube();
    void cube_based_estimate();
    void iterate_through_all_cubes(int current_d);
};

#endif //FINAL_CS206_CUBEBASED_H

#pragma once
#include <Eigen/Dense>//<eigen/3.4.0_1/include/eigen3/Eigen/Dense>

using namespace Eigen;

class Grid_Point{
    public:
        Grid_Point();
        ~Grid_Point();
        void Init();

    public:
    double mass;
    double active;
    bool fix_tag;
    Vector3d velocity;
    Vector3d velocity_col;
    Vector3d velocity_fri;
    Vector3d momentum;
    Vector3d force;
    Vector3d g;
};
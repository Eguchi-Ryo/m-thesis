#pragma once
#include <Eigen/Dense>//<Eigen/3.4.0_1/include/eigen3/Eigen/Dense>
#include <vector>
#include "Parameters.h"

using namespace std;
using namespace Eigen;
enum MaterialType {PANC, STOM, DUO, LKID};

class Particle {
    public:
    Particle();
    ~Particle();

    Particle(Vector3d Position, Vector3d Velocity, double Mass, double Density, double Dcmlabel, MaterialType Tag, unsigned short Dcmval, int Z_index, bool Force_tag);

    //関数
 
    void Compute_Stress(); //応力を計算

    void Apply_Damping();
    Vector3d Forced_Stress();
    Vector3d Damping_Force();

    //変数
    double mass;
    double volume;
    double init_volume;
    double density;
    double init_density;
    double J;



    MaterialType tag;
    int force_tag;
    bool fix_tag;


    Vector3d position;
    Vector3d pseudo_position;
    Vector3i int_position;
    Vector3d init_position;
    Vector3d target_position;
    Vector3d velocity;
    Vector3d momentum;
    Vector3d weight_gradient[4][4][4];
    Vector3d bef_position;
    Vector3d integral;

    Matrix3d F_e; //変形勾配の弾性成分
    Matrix3d L; //速度勾配
    Matrix3d C; //アフィン速度行列
    Matrix3d Energy_Derivative; //応力
    vector<vector<vector<double > > > weight;

    //形状関数
    Matrix4d type_x;
    Matrix4d type_y;
    Matrix4d type_z;

    double dcmval;
    double dcmlabel;
    int z_index;

    double prop = 1e10;//1e8
    double diff = prop * 0.8;//2.5e6 
    double inte = prop * 0.1;//5e6 

    Vector3d error_buf_pid;

    


};
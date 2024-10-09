#include "../Header/Particle.h"
#include "../Header/Parameters.h"
#include <cstdlib>
#include <iostream>
#include "../Header/Grid.h"
#include <algorithm>
using namespace Eigen;

#define MAX_PID_FORCE 1e11
Particle::Particle(){

}

Particle::~Particle(){

}

Particle::Particle(Vector3d Position, Vector3d Velocity, double Mass, double Density, double Dcmlabel, MaterialType Tag, unsigned short Dcmval, int Z_index, bool Force_tag){
    position = Position;
    velocity = Velocity;
    mass = Mass;
    momentum = mass * velocity;
    tag = Tag;
    density = Density;
    init_density = Density;
    volume = mass / density;
    init_volume = mass / density;
    dcmlabel = Dcmlabel;
    Energy_Derivative = Matrix3d::Zero();
    J = 1; //　変形勾配全体の行列式
    C = Matrix3d::Zero();
    F_e = Matrix3d::Identity();

    //cout << "finish" << endl;
    weight.resize(4);
    for(int i = 0; i < 4; i ++){
        weight[i].resize(4);
            for(int j = 0; j < 4; j++){
                weight[i][j].resize(4);
                for(int k = 0; k < 4; k++)weight[i][j][k] = 0;
            }
    }
    fix_tag = false;
    force_tag = Force_tag;
    init_position = position;
    target_position = position;
    bef_position = position;
    integral = Vector3d::Zero();
    dcmval = Dcmval;
    z_index = Z_index;
    error_buf_pid = Vector3d::Zero();

}

 
void Particle::Compute_Stress(){

    JacobiSVD<MatrixXd>svd(F_e,ComputeThinU | ComputeThinV);
    Matrix3d U = svd.matrixU(); //左直交行列
    Matrix3d V = svd.matrixV(); //右直交行列
    Matrix3d S = U.transpose() * F_e * V; //ひずみテンソルの回転不変量
    double det = S.determinant();

    if(tag == PANC){

        Energy_Derivative = 2 * PANC_MIU * (F_e - U * V.transpose()) * F_e.transpose() + PANC_LAMBDA * (det - 1) * det * Matrix3d::Identity();
    
    }else if(tag == STOM){

        Energy_Derivative = 2 * STOM_MIU * (F_e - U * V.transpose()) * F_e.transpose() + STOM_LAMBDA * (det - 1) * det * Matrix3d::Identity();
       

    }else if(tag == DUO){

        Energy_Derivative = 2 * DUO_MIU * (F_e - U * V.transpose()) * F_e.transpose() + DUO_LAMBDA * (det - 1) * det * Matrix3d::Identity();
       
    }else if(tag == LKID){

        Energy_Derivative = 2 * LKID_MIU * (F_e - U * V.transpose()) * F_e.transpose() + LKID_LAMBDA * (det - 1) * det * Matrix3d::Identity();
       
    }
    else{
        //何もしない
    }
}

Vector3d Particle::Damping_Force(){
    double damping_force = 0;
    Vector3d out;
    
    if(force_tag != 1){
        
        //各パラメータの小さい方の1/10を係数とする
        if(tag == PANC){
            //damping_force = min(PANC_LAMBDA,PANC_MIU) * 0.1;
            damping_force = 2 * sqrt(3.5e3);
        }else if(tag == STOM){
            //damping_force = min(STOM_LAMBDA,STOM_MIU) * 0.1;
            damping_force = 2 * sqrt(500e3);
        }else if(tag == DUO){
            //damping_force = min(DUO_LAMBDA,DUO_MIU) * 0.1;
            damping_force = 2 * sqrt(1e7);
        }else if(tag == LKID){
            //damping_force = min(LKID_LAMBDA,LKID_MIU) * 0.1;
            damping_force = 2 * sqrt(24e3);
        }
        
        //damping_force = 1000;
        //cout << damping_force << endl;
        out << damping_force / DT * (position[0] - bef_position[0]),
               damping_force / DT * (position[1] - bef_position[1]),
               damping_force / DT * (position[2] - bef_position[2]);
    }else out = Vector3d::Zero();

    return -out;
}

//強制変位を加えた際の力
Vector3d Particle::Forced_Stress(){
    Vector3d out;

    Vector3d now_error = position - target_position;

    if(force_tag == 1){
        integral += (error_buf_pid + now_error) / 2.0 * DT;
        
        out <<  prop * now_error[0] + diff * (now_error[0] - error_buf_pid[0]) / DT + inte * integral[0],
                prop * now_error[1] + diff * (now_error[1] - error_buf_pid[1]) / DT + inte * integral[1], 
                prop * now_error[2] + diff * (now_error[2] - error_buf_pid[2]) / DT + inte * integral[2];
        
        //cout << out << endl << endl;
        /*
        integral[0] += (position[0] + bef_position[0] - 2 * target_position[0]) / 2.0 * DT;
        integral[1] += (position[1] + bef_position[1] - 2 * target_position[1]) / 2.0 * DT;
        integral[2] += (position[2] + bef_position[2] - 2 * target_position[2]) / 2.0 * DT;
        out <<  prop * (position[0] - target_position[0]) + diff / DT * (position[0] - bef_position[0]) + inte * integral[0],
                prop * (position[1] - target_position[1]) + diff / DT * (position[1] - bef_position[1]) + inte * integral[1], 
                prop * (position[2] - target_position[2]) + diff / DT * (position[2] - bef_position[2]) + inte * integral[2];
        */
       error_buf_pid = now_error;   
    

        
    }else out = Vector3d::Zero();
    
    // 大きさの最大値に範囲を抑える
    for(int i = 0; i < 3; i++){
        if(out(i) >= MAX_PID_FORCE){
            out(i) = MAX_PID_FORCE;
        }else if(out(i) <= -1 * MAX_PID_FORCE){
            out(i) = -1 * MAX_PID_FORCE;
        }
    }
    
    return -out;

}



//ダンピング処理
void Particle::Apply_Damping(){
    Matrix3d skew, symm;
    skew = C - C.transpose();
    skew *= 0.5;
    symm = C - skew;
    C = skew + (1-DAMPRATIO)*symm;
}
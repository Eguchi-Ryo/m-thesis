#pragma once
#include "Grid_Point.h"
#include "Particle_Cloud.h"
#include <vector>
#include <Eigen/Dense>//<eigen/3.4.0_1/include/eigen3/Eigen/Dense>
#include <Eigen/Geometry>//<eigen/3.4.0_1/include/eigen3/Eigen/Geometry>

using namespace Eigen;
using namespace std;
extern int Duodenum_num;
extern int Pancreas_num;
extern int Stomach_num;
extern int Left_Kidney_num;

extern bool pid_output_flag;
extern int mode;
extern bool unconverged_flag;
extern bool divergence_flag;

extern string th_csv;
extern double TH;

extern string affine3D;
extern string df_csv;
extern string cs_list;
extern string df_ratio_list;

extern string deformationField_3D;

extern string trans_fp; // 変換パラメータファイルへのパス
extern int cs; //力を加える断面

/*複数方向の断面用（直接画素値が存在するものを指定する）*/
extern int cs_x;
extern int cs_y;
extern int cs_z;
extern ofstream logout;
extern string target_csv; //ターゲットの目標位置のファイルパス

extern double img_ratio[3]; // 画像の倍率

extern int image_size[3];

extern double pixel_spacing[3];

extern double DT;

class Grid{
    public:
        Grid();
        ~Grid();
        Grid(Vector3d Start, Vector3d End, Vector3i Size, Particle_Cloud* Object);
        //double N(double x);
        //double dN(double x); //Nの微分
        double N(double x,int type);
        double dN(double x,int type); //Nの微分

        //通常のアルゴリズム
        void Init_Grid();

        void Update_From_Material();
        void Update_Grid();
        void Compute_Weight_Gradient();
        void Reset_Force();
        void Update_Force();
        void Compute_Velocity();
        void Determine_Collision();
        void Update_Particle_Cloud();
        void Update_Particle_Position();
        void Update_Deformation_Gradient();


        

        //タイムステップ計算用の関数
        vector<double> Calc_C();
        void Calc_DT();

        //断面制御用の関数
        void Make_Target_Position();
        void Output_Elastix_Position();
        void Make_Target_3D_Position();
        void Calc_Error();
        void Calc_Threshold();

        void Make_Target_multi_cs_Position();



    public:
        vector<vector<vector<Grid_Point> > > Map;
        Vector3i size; //格子点数
        Vector3d start, end; //格子点の始点、終点
        Vector3d cellsize; //各格子点が保持する領域のサイズ
        Particle_Cloud* Material; //対象の情報
        int constraint; //端から何個の点が境界線であるかどうかを決める変数
        int step_g; //格子点のステップ数
        double error; //目標位置との誤差
        double error_buf;

        double early_stop;
        double early_stop2;
        double init_error;
        double error_dd;
        double error_dd_buf;
        bool error_flag;
        double error_flag_num;

        vector<double> sound;
        vector<double> max_vel;

        double pid_dump;
        Vector3d Max_PID_Force;
        Vector3d Min_PID_Force;
        Vector3d Max_Dumping_Force;






};


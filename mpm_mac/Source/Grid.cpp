#include "../Header/Grid.h"
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include <omp.h>///</opt/homebrew/opt/libomp/include/omp.h>
#include <opencv2/opencv.hpp>
#include <iomanip>
#include <array>

#define ERROR_EPS 1e-6
#define ERROR_STEPS 100
#define MAX_STEPS 2000
#define EARLY_STOP_NUM 100
#define REL_ERROR_EPS 0.0015

int Stomach_num = 0;
int Pancreas_num = 0;
int Duodenum_num = 0;
int Left_Kidney_num = 0;
int mode = 2;
bool pid_output_flag = false;
bool unconverged_flag = false;
bool divergence_flag = false;
string trans_fp = "";
double TH = INIT_THRESHOLD;
string th_csv = "";
ofstream logout;
int cs = 0;
int cs_x = -1;
int cs_y = -1;
int cs_z = -1;
vector<vector<double> > target_info;
string target_csv = "";
double img_ratio[3] = {0,0,0};
double DT = 0;
//double res_error = 0;
Grid::Grid(){

}

Grid::~Grid(){
    delete Material;
}

Grid::Grid(Vector3d Start, Vector3d End, Vector3i Size, Particle_Cloud* Object){
    start = Start;
    end = End;
    size = Size;
	
    constraint = 3; //境界条件を決める際に使う。端から何個のセルが境界かどうかを表す
    cellsize[0] = (end[0] - start[0]) / (double)(size[0] - 1);
    cellsize[1] = (end[1] - start[1]) / (double)(size[1] - 1);
	cellsize[2] = (end[2] - start[2]) / (double)(size[2] - 1);

    Map.resize(Size[0]);
    for(int i = 0; i < Size[0]; i++){
        Map[i].resize(Size[1]);
		for(int j = 0; j < Size[1]; j++){
			Map[i][j].resize(Size[2]);
		}
    }
    Material = Object;
    step_g = 0;
	error = 0;
	error_buf = 0;
	init_error = 0;
	error_dd = 0;
	error_dd_buf = 0;
	error_flag_num = 0;
	error_flag = false;

	//mode1:データ作成処理
	//mode2:推定処理
	cout << "TH : " << TH << endl;
	cout << "mode = " << mode << endl;
	if (mode == 1) {
		Make_Target_3D_Position();
	}
	if (mode == 2) {
		Make_Target_Position();
	}
	if (mode == 4){
		Make_Target_multi_cs_Position();
	}


	
	//Output_Elastix_Position();

	sound.resize(3);
	sound = Calc_C();
	max_vel.resize(3,0);
	Calc_DT();
	pid_dump = 1;
	Max_PID_Force = Vector3d::Zero();
	Min_PID_Force = Vector3d::Zero();
	Min_PID_Force << 1e20, 1e20, 1e20;
	Max_Dumping_Force = Vector3d::Zero();
	early_stop = 0;
	early_stop2 = 0;
	

}

/***********************************************************************************************************************/
/***********************************************断面から制御をするための関数群***********************************************/
/***********************************************************************************************************************/

//画像を変形させ、目的の座標を決定する
// 座標変換関数
//arrayの１ボクセルに対する処理
std::array<double, 3> adjustCoordinateSystem(const std::array<double, 3>& displacement) {
    // 左手系から右手系への変換
	return { -displacement[0], displacement[1], displacement[2] };
    //return { displacement[0], displacement[1], displacement[2] };  //patient, (A,S,C)=(96, 200, 180)の時完全一致
}

void Grid::Make_Target_Position(){
    VectorXd vec = VectorXd::Zero(4);
    VectorXd new_vec = VectorXd::Zero(4);
	VectorXd tar_SimSpace = VectorXd::Zero(4);
	string DF_RawFile = deformationField_3D;
    // データのサイズ情報
    int width = image_size[0];
    int height = image_size[1];
    int depth = image_size[2];

	/*
    std::vector<std::vector<std::vector<std::array<double, 3>>>> deformationField(
        depth, std::vector<std::vector<std::array<double, 3>>>(
            height, std::vector<std::array<double, 3>>(width)));
	*/
    std::vector<std::vector<std::vector<std::array<double, 3>>>> deformationField(
        width, std::vector<std::vector<std::array<double, 3>>>(
            height, std::vector<std::array<double, 3>>(depth)));

    // RAWバイナリファイルを読み込む
    std::ifstream file(DF_RawFile, std::ios::binary);
    if (!file) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // ファイルからデータを3次元配列に読み込む
	//cppはrow-major orderであることに留意
    for (int z = 0; z < depth; ++z) {
	//for (int z = depth - 1; z >= 0; --z) {
		//for (int y = height - 1; y >= 0; --y) { 
		for (int y = 0; y < height; ++y) {
		//for (int x = 0; x < width; ++x) { 
			for (int x =0; x < width; ++x) {
				std::array<double, 3> displacement;
	
				file.read(reinterpret_cast<char*>(&displacement), sizeof(double) * 3);
				//deformationField[depth -1 -z][height -1 -y][x] = adjustCoordinateSystem(displacement);  
				deformationField[x][height -1 -y][z] = adjustCoordinateSystem(displacement); 
				
			}
		}
	}



    auto displacement2 = deformationField[depth -1 -66][height -1 -244][283];
    std::cout << "Displacement at [66, 244, 283], (dx,dy,dz)= (" 
              << displacement2[0] << ", " 
              << displacement2[1] << ", " 
              << displacement2[2] << ")" << std::endl;

    auto displacement3 = deformationField[232][0][0];
    std::cout << "Displacement at [233, 0, 0], (dx,dy,dz)= (" 
              << displacement3[0] << ", " 
              << displacement3[1] << ", " 
              << displacement3[2] << ")" << std::endl;
	
	
	std::cout << "displacement size " << deformationField.size() << endl;

    file.close();
	cout << "complete reading deformationField file" << endl;

    for(int p = 0; p < Material->points.size(); p++){
        Particle& P = Material->points[p];
		P.force_tag = 1;

		int pos_x = P.position[0] / img_ratio[0];
		int pos_y = P.position[1] / img_ratio[1];
		int pos_z = P.position[2] / img_ratio[2];

        vec(0) = P.position[0] / img_ratio[0] * pixel_spacing[0];  
        vec(1) = P.position[1] / img_ratio[1] * pixel_spacing[1]; 
        vec(2) = P.position[2] / img_ratio[2] * pixel_spacing[2];
        vec(3) = 1;

		//std::array<double, 3> disp = deformationField[pos_z][pos_y][pos_x];   
		std::array<double, 3> disp = deformationField[pos_x][pos_y][pos_z];  

		float dx_sim = disp[0]*img_ratio[0];
		float dy_sim = disp[1]*img_ratio[1];
		float dz_sim = disp[2]*img_ratio[2];

		tar_SimSpace(0) = P.position[0] + dx_sim;
		tar_SimSpace(1) = P.position[1] + dy_sim;
		tar_SimSpace(2) = P.position[2] + dz_sim; 
		
		new_vec(0) = vec(0) + dx_sim / img_ratio[0] * pixel_spacing[0];
		new_vec(1) = vec(1) + dy_sim / img_ratio[1] * pixel_spacing[1]; 
		new_vec(2) = vec(2) + dz_sim / img_ratio[2] * pixel_spacing[2];

        P.target_position[0] = new_vec(0) * img_ratio[0] / pixel_spacing[0];
        P.target_position[1] = new_vec(1) * img_ratio[1] / pixel_spacing[1]; 
        P.target_position[2] = new_vec(2) * img_ratio[2] / pixel_spacing[2];
		
	}
	cout << "done particle loop" << endl;
}

//3次元変位
//csvから各質点の目標位置を計算する
void Grid::Make_Target_3D_Position(){
	string fname = target_csv;
	string trans_fname = trans_fp;
	string buf, buf2, param;
	int counter = 0;
	double buf_d;


	vector<double> array;
	ifstream read_file(target_csv);
	//csvファイルから目標位置を読み込み
	while(getline(read_file,buf)){
		istringstream input(buf);
		counter = 0;
		while(getline(input,buf2,',')){
			buf_d = stod(buf2);
			array.push_back(buf_d);
			//cout << buf2 << " ";
		}
		//cout << endl;
		target_info.push_back(array);

		array.clear();
	}



	cout << "calc target position " << endl;
	for(int p = 0; p < Material->points.size(); p++){
		Particle& P = Material->points[p];

		P.force_tag = 1;
		//cout << target_info.size() << " " << Material->points.size() << endl;
		if(target_info[p][6] != P.dcmlabel){
			cout << "doesn't same !" << endl;
			cout << target_info[p][5] << " " << target_info[p][6] << " " << P.dcmval << " " << P.dcmlabel << endl;
			exit(1);
		}
		if(target_info[p][5] != P.dcmval) {
			cout << "doesn't same !" << endl;
			cout << target_info[p][5] << " " << target_info[p][6] << " " << P.dcmval << " " << P.dcmlabel << endl;
			exit(1);
		}
		P.target_position[0] = target_info[p][1] * img_ratio[0];
		P.target_position[1] =target_info[p][2] * img_ratio[1];
		P.target_position[2] = target_info[p][3] * img_ratio[2];
          
	}
	cout << "image ratio in Grid: img_raio[0]" << img_ratio[0] << "img_ratio[1]" << img_ratio[1] << "img_ratio[2]" << img_ratio[2] << endl;

}

//複数方向のcsを使用する場合
//csvから各質点の目標位置を計算する
void Grid::Make_Target_multi_cs_Position(){
	string fname = target_csv;
	string trans_fname = trans_fp;
	string buf, buf2, param;
	int counter = 0;
	double buf_d;


	vector<double> array;
	ifstream read_file(target_csv);
	//csvファイルから目標位置を読み込み
	while(getline(read_file,buf)){
		istringstream input(buf);
		counter = 0;
		while(getline(input,buf2,',')){
			buf_d = stod(buf2);
			array.push_back(buf_d);
			//cout << buf2 << " ";
		}
		//cout << endl;
		target_info.push_back(array);

		array.clear();
	}

	int particle_sum = 0;
	cout << "calc cs target position " << endl;

	for(int p = 0; p < Material->points.size(); p++){
		Particle& P = Material->points[p];
		//if(P.z_index == cs || int(P.position[0]/img_ratio[0]) == x_thr || int(P.position[1] / img_ratio[1]) == y_thr){
		//if(int(P.position(0)) == x_thr || int(P.position(1)) == y_thr){
		if(P.force_tag){
			//P.force_tag = 1;
			//cout << P.position << endl;
			particle_sum++;

			if(target_info[p][6] != P.dcmlabel){
				cout << "doesn't same !" << endl;
				cout << target_info[p][5] << " " << target_info[p][6] << " " << P.dcmval << " " << P.dcmlabel << endl;
				exit(1);
			}
			if(target_info[p][5] != P.dcmval) {
				cout << "doesn't same !" << endl;
				cout << target_info[p][5] << " " << target_info[p][6] << " " << P.dcmval << " " << P.dcmlabel << endl;
				exit(1);
			}
			P.target_position[0] = target_info[p][1] * img_ratio[0];
			P.target_position[1] =target_info[p][2] * img_ratio[1];
			P.target_position[2] = target_info[p][3] * img_ratio[2];
		}

	}
	cout << "particle sum : " << particle_sum << endl;
	cout << "target_size " << target_info.size() << endl;

}

//制御断面との誤差を計算。閾値を下回った場合にフラグを立てる
void Grid::Calc_Error(){
	//double TH = 0.005;
	double disp_TH = 0.5;
	error_buf = error;
	error_dd_buf = error_dd;
	error = 0;
	int num = 0;
	for(int p = 0; p < Material->points.size(); p++){
		if(Material->points[p].force_tag && Material->points[p].dcmlabel == 1){
			//error += (Material->points[p].position - Material->points[p].target_position).norm();
			error += (Material->points[p].position - Material->points[p].target_position).norm() /(Material->points[p].target_position).norm();
			num++;
		}
	}
	//cout << error << endl;
	error /= (double)num;
	//cout << "mother: " << (double)num << endl;
	//res_error /= (double)num;

	cout << "error = " << error << endl;
	logout << "error = " << error << endl;
	/*
	if(step_g == 1){
		init_error = error;
		cout << "init_error : " << init_error << endl;
	}
	error_dd = error - error_buf;
	cout << "error_dd = " << error_dd << endl;
	if (step_g > 100 && error_flag_num == 0 && abs(error - error_buf) <= ERROR_EPS){
		error_flag_num++;
	}else if(error_flag_num != 0 && abs(error - error_buf) <= ERROR_EPS){
		error_flag_num++;
	}else{
		error_flag_num = 0;
	}

	if (step_g > 100 && abs(error - error_buf) <= ERROR_EPS && error < error_buf){
		early_stop++;
	}else{
		early_stop = 0;
	}*/
	if (error < REL_ERROR_EPS){
		early_stop++;
	}else{
		early_stop = 0;
	}

	//最高位とその次の桁を求める
	string ebuf_s = to_string(error_buf);
	string e_s = to_string(error);
	ebuf_s = ebuf_s.erase(6);
	e_s = e_s.erase(6);
	//cout << ebuf_s << " " << e_s << endl; //元々出してた
	if(ebuf_s==e_s){
		early_stop2++;
	}else{
		early_stop2 = 0;
	}
	 //元々出していた
	cout << "error_diff : " << abs(error - error_buf) << endl;
	cout << "early_stop num : " << early_stop << endl;
	cout << "early_stop2 num : " << early_stop2 << endl;
	//cout << "error_flag_num : " << error_flag_num << endl;
	cout << "error_flag : " << error_flag << endl;

	/*
	if((error_flag && error_dd > 0 && error_dd_buf < 0) || early_stop >= EARLY_STOP_NUM){
		pid_output_flag = true;
	}else if(error_flag_num == ERROR_STEPS){// || (error_buf - error < 0 && error < TH*2)) {
		error_flag = true;
		//pid_output_flag = true;
		//getchar();
	}
	*/
	if (early_stop == EARLY_STOP_NUM || early_stop2 == EARLY_STOP_NUM){
		pid_output_flag = true;
	}else if(step_g >= MAX_STEPS){ 
	//THを下回らなかった場合はdataフォルダに番号と条件を出力する。
	//最初の数ステップでは多少増加する場合があるので、ある程度step数が増えてからということにする。
		cout << "the error does not converge well" << endl;
		logout << "the error does not converge well" << endl;
		pid_output_flag = true;
		unconverged_flag = true;
		//exit(0);
	}else if(error > disp_TH || error_buf == error){
		cout << "simulation divergence" << endl;
		logout << "simulation divergence" << endl;
		divergence_flag = true;
		pid_output_flag = true;
		//getchar();
		//exit(1);
	}


	
	if(step_g > 10000){ 
		pid_output_flag = true;
	}
	
	
}


/***********************************************************************************************************************/
/***************************************************CFL条件のための計算****************************************************/
/***********************************************************************************************************************/

vector<double> Grid::Calc_C(){
	//胃、十二指腸の大きい方のcを返す。二つあるのは格子を二つ定義した時用。膵臓と腎臓は小さいので除外
	double stom_young = STOM_MIU*(3*STOM_MIU + 2 * STOM_LAMBDA) / (STOM_MIU + STOM_LAMBDA);
	double stom_poisson = STOM_LAMBDA / (2 * (STOM_LAMBDA + STOM_MIU));
	double duo_young = DUO_MIU*(3*DUO_MIU + 2 * DUO_LAMBDA) / (DUO_MIU + DUO_LAMBDA);
	double duo_poisson = DUO_LAMBDA / (2 * (DUO_LAMBDA + DUO_MIU));

	//密度を変えたら変更する必要あり
	double c1 = sqrt(stom_young * (1 - stom_poisson) / ((1 + stom_poisson) * (1 - 2 * stom_poisson) * 1));
	double c2 = sqrt(duo_young * (1 - duo_poisson) / ((1 + duo_poisson) * (1 - 2 * duo_poisson) * 1));
	
	double c_max = max(c1,c2);
	vector<double> c(2);
	c[0] = c_max ;
	c[1] = c_max ;
	return c ;

}

//CFL条件を元にタイムステップを計算
void Grid::Calc_DT(){
	double alpha = 0.4; //0.5

	double t_x1 = 1.0 / (size[0] * max(sound[0],max_vel[0]));
	double t_y1 = 1.0 / (size[1] * max(sound[0],max_vel[1]));
	double t_z1 = 1.0 / (size[2] * max(sound[0],max_vel[2]));

	double t_x2 = 0;
	double t_y2 = 0;
	double t_z2 = 0;

	double t_x = max(t_x1,t_x2);
	double t_y = max(t_y1,t_y2);
	double t_z = max(t_z1,t_z2);
	//cout << "tx : " << t_x << " ty : " << t_y << " tz : " << t_z << endl;  //元々出していた
	//PID的に力を加えるため、普通のCFL条件では計算が発散してしまったため小さくした
	//DT = alpha * max(max(t_x, t_y),t_z) * 0.1;
	
	DT = alpha * min(min(t_x, t_y),t_z);

}


/************************************************************************************************************************/
/*****************************************************MPMの為の関数群******************************************************/
/************************************************************************************************************************/


void Grid::Reset_Force(){
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			for(int k = 0; k < size[2]; k++){
				Map[i][j][k].force = Vector3d::Zero();
			}
		}
	}
}

/*
//枠との関係を考えない場合
double Grid::N(double x){
	double N_x;
	x = fabs(x);
	if (x < 1)
		N_x = x * x * (x / 2 - 1) + 2 / 3.0;
	else if (x < 2)
		N_x = x * (x * (-x / 6 + 1) - 2) + 4 / 3.0;
	else
		N_x = 0.0;

	if(N_x < BSPLINE_EPSILON) return 0;
	return N_x;
}
*/

//枠との関係で補間関数を変える場合
double Grid::N(double x, int type){
	double N_x;
	if(type == 1){
		x = fabs(x);
		if (x < 1)
			N_x = x * ( x * x / 6 - 1) + 1.0;
		else if (x < 2)
			N_x = x * (x * (-x / 6 + 1) - 2) + 4 / 3.0;
		else
			N_x = 0.0;		
	}else if(type == 2){
		if ( -1 <= x && x <= 0)
			N_x = -1 / 3.0 * x * x * x - x * x + 2 / 3.0;
		else if (0 <= x && x <= 1)
			N_x = 1 / 2.0 * x * x * x - x * x + 2 / 3.0;
		else if (1 <= x && x <= 2)
			N_x = -1/6.0 * x * x * x + x * x - 2 * x + 4 / 3.0;
		else 
			N_x = 0.0;
	}else if(type == 3){
		x = fabs(x);
		if (x < 1)
			N_x = x * x * (x / 2 - 1) + 2 / 3.0;
		else if (x < 2)
			N_x = x * (x * (-x / 6 + 1) - 2) + 4 / 3.0;
		else
			N_x = 0.0;		
	}else if(type == 4){
		if (-2 <= x && x <= -1)
			N_x = 1/6.0 * x * x * x + x * x + 2 * x + 4 / 3.0;
		else if (-1 <= x && x <= 0)
			N_x = - 1 / 2.0 * x * x * x - x * x + 2 / 3.0;
		else if (0<= x && x <= 1)
			N_x = 1 / 3.0 * x * x * x - x * x + 2 / 3.0;
		else 
			N_x = 0.0;
	}
	
	if(N_x < BSPLINE_EPSILON) return 0;
	return N_x;
}

/*
//枠との関係を考えない場合
double Grid::dN(double x){

	double abs_x = fabs(x);
	if (abs_x < 1)
		return 1.5*x*abs_x - 2 * x;
	else if (x < 2)
		return -x * abs_x / 2 + 2 * x - 2 * x / abs_x;
	else return 0;
}
*/

//枠との関係で補間関数を変える場合
double Grid::dN(double x,int type){

	if(type == 1){
		double abs_x = fabs(x);
		if (abs_x < 1)
			return 0.5*x*abs_x - abs_x;
		else if (x < 2)
			return -x * abs_x / 2 + 2 * x - 2 * x / abs_x;
		else return 0;		
	}else if (type == 2){
		if ( -1 <= x && x <= 0)
			return - x * x - 2 * x;
		else if (0 <= x && x <= 1)
			return 3 / 2.0 * x * x - 2 * x;
		else if (1 <= x && x <= 2)
			return - 0.5 * x * x + 2 * x - 2;
		else 
			return 0.0;		
	}else if (type == 3){
		double abs_x = fabs(x);
		if (abs_x < 1)
			return 1.5*x*abs_x - 2 * x;
		else if (x < 2)
			return -x * abs_x / 2 + 2 * x - 2 * x / abs_x;
		else return 0;
	}else if (type == 4){
		if (-2 <= x && x <= -1)
			return x * x + 2 * x + 2;
		else if ( x <= 0)
			return - 3.0 / 2.0 * x * x - 2 * x;
		else if ( x <= 1)
			return x * x - 2 * x;
		else 
			return 0.0;
	}else{
		return 0;
	}

}


void Grid::Update_From_Material(){
	//cout << "b" << endl;
	#pragma omp parallel for
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			for(int k = 0; k < size[2]; k++){
				Map[i][j][k].mass = 0;
				Map[i][j][k].velocity = Vector3d::Zero();
				Map[i][j][k].active = false;
				Map[i][j][k].momentum = Vector3d::Zero();
				Map[i][j][k].fix_tag = false;
			}

		}
	}
	for (int p = 0; p < Material->points.size(); p++)
	{
		Particle& P = Material->points[p];

        //格子点の通し番号に変換→整数にする
		double x = P.position[0], y = P.position[1], z = P.position[2];
		x = (x - start[0]) / cellsize[0], y = (y - start[1]) / cellsize[1], z = (z - start[2]) / cellsize[2];
		P.pseudo_position = Vector3d(x, y, z);
		int u = (int)floor(x), v = (int)floor(y), r = (int)floor(z);
		P.int_position = Vector3i(u, v, r);

		#pragma omp parallel for
		for(int i = 0; i < 4; i++){
			for (int j = 0; j < 4; j++){
				for(int k = 0; k < 4; k++){
					if (i + u -1 >= 0 && (i + u - 1) < size[0] && j + v - 1 >= 0 && (j + v - 1) < size[1] && k + r - 1 >= 0 && (k + r - 1) < size[2]){
						//境界条件
						if(P.fix_tag) Map[i + u - 1][j + v - 1][k + r - 1].fix_tag = true;
					}
				}
			}
		}

		#pragma omp parallel for
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				for(int k = 0; k < 4; k++){
					if (i + u -1 >= 0 && (i + u - 1) < size[0] && j + v - 1 >= 0 && (j + v - 1) < size[1] && k + r - 1 >= 0 && (k + r - 1) < size[2]) {
						//固定を確認
						for(int b = -1; b <= 1; b++){
							if(i + u - 1 + b >= 0 && i + u - 1 + b < size[0]) {
								if(Map[i + u -1 + b][j + v -1][k + r - 1].fix_tag && b == 0) P.type_x(i,j) =  1;
								else if(Map[i + u -1 + b][j + v -1][k + r - 1].fix_tag && b == -1) P.type_x(i,j) = 2;
								else if(Map[i + u -1 + b][j + v -1][k + r - 1].fix_tag && b == 1) P.type_x(i,j) = 4; 
								else P.type_x(i,j) = 3;
							}

							if(j + v - 1 + b >= 0 && j + v - 1 + b < size[1]) {
								if(Map[i + u -1][j + v -1 + b][k + r - 1].fix_tag && b == 0) P.type_y(i,j) =  1;
								else if(Map[i + u -1][j + v -1 + b][k + r - 1].fix_tag && b == -1) P.type_y(i,j) = 2;
								else if(Map[i + u -1][j + v -1 + b][k + r - 1].fix_tag && b == 1) P.type_y(i,j) = 4; 
								else P.type_y(i,j) = 3;
							}
							if(k + r - 1 + b >= 0 && k + r - 1 + b < size[2]) {
								if(Map[i + u -1][j + v -1][k + r - 1 + b].fix_tag && b == 0) P.type_z(i,j) =  1;
								else if(Map[i + u -1][j + v -1][k + r - 1 + b].fix_tag && b == -1) P.type_z(i,j) = 2;
								else if(Map[i + u -1][j + v -1][k + r - 1 + b].fix_tag && b == 1) P.type_z(i,j) = 4; 
								else P.type_z(i,j) = 3;
							}
						}
						P.weight[i][j][k] = N(i + u - x - 1,P.type_x(i,j)) * N(j + v - y - 1,P.type_y(i,j))* N(k + r - z - 1,P.type_z(i,j));
						if (P.weight[i][j][k] > BSPLINE_EPSILON)
							Map[i + u - 1][j + v - 1][k + r - 1].active = true;
							#pragma omp critical
							Map[i + u - 1][j + v - 1][k + r - 1].mass += P.weight[i][j][k] * P.mass;
							//cout << Map[i + u - 1][j + v - 1][k + r - 1].mass << endl;
							Eigen::Vector3d position = Eigen::Vector3d(start[0] + cellsize[0] * (i + u - 1), start[1] + cellsize[1] * (j + v - 1), start[2] + cellsize[2] * (k + r - 1));
							#pragma omp critical
							//運動量を求める
							Map[i + u - 1][j + v - 1][k + r - 1].momentum += P.weight[i][j][k] * P.mass*(P.velocity + 3.0 / cellsize[0] / cellsize[1] * P.C*(position - P.position));

					}
				}
			}
		}
	}
	#pragma omp parallel for
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			for(int k = 0; k < size[2]; k++){
				if (Map[i][j][k].active && Map[i][j][k].fix_tag == false && Map[i][j][k].mass > MASS_EPSILON){
					Map[i][j][k].velocity = Map[i][j][k].momentum / Map[i][j][k].mass;//運動量を質量で割ることで速度を求める
					//cout << "vel " << Map[i][j][k].velocity << endl;
				}else {
					Map[i][j][k].velocity = Vector3d(0, 0, 0);
					Map[i][j][k].momentum = Vector3d(0,0,0);
				}
			}

		}
	}
	
}

//重み関数の空間微分を計算
//形状関数はN(x_vector)=N(1/h(x-ih))*N(1/h(y-ih))で表される為グリッド間距離hで式を割っている
void Grid::Compute_Weight_Gradient(){
	
	for (int p = 0; p < Material->points.size(); p++)
	{
		Particle& P = Material->points[p];
		double x = P.pseudo_position[0], y = P.pseudo_position[1], z = P.pseudo_position[2];
		int u = P.int_position[0], v = P.int_position[1], r = P.int_position[2];
		#pragma omp parallel for
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				for(int k = 0; k < 4; k++){
					P.weight_gradient[i][j][k] << dN(x - u - i + 1,P.type_x(i,j)) * N(y - v - j + 1,P.type_y(i,j))*N(z - r - k + 1,P.type_z(i,j)) / cellsize[0],
					N(x - u - i + 1,P.type_x(i,j)) * dN(y - v - j + 1,P.type_y(i,j)) * N(z - r - k + 1,P.type_z(i,j)) / cellsize[1], 
					N(x - u - i + 1,P.type_x(i,j)) * N(y - v - j + 1,P.type_y(i,j))*dN(z - r - k + 1,P.type_z(i,j)) / cellsize[2];
				}

			}
		}
	}
}

void Grid::Update_Force(){
	Reset_Force();
	for (int p = 0; p < Material->points.size(); p++)
	{
		Particle& P = Material->points[p];
		double x = P.pseudo_position[0], y = P.pseudo_position[1], z = P.pseudo_position[2];
		int u = P.int_position[0], v = P.int_position[1], r = P.int_position[2];
		#pragma omp parallel for
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				for(int k = 0; k < 4; k++){
					if (u + i - 1 >= 0 && u + i - 1 < size[0] && v + j - 1 >= 0 && v + j - 1 < size[1] && k + r - 1 >= 0 && (k + r - 1) < size[2]) {
						//Map[u + i - 1][v + j - 1][r + k - 1].force = Map[u + i - 1][v + j - 1][r + k - 1].force - P.volume * P.Energy_Derivative * P.weight_gradient[i][j][k] + P.weight[i][j][k] * P.Forced_Stress() * P.volume / cellsize[2];
						//ダンピング力を加える場合
						Vector3d PID_force = P.Forced_Stress();
						Vector3d Dumping_force = P.Damping_Force();
						/*
						if( p == 0) {
							cout << P.position - P.target_position << endl << endl;
							cout << P.error_buf << endl;
						}
						*/
						Map[u + i - 1][v + j - 1][r + k - 1].force = Map[u + i - 1][v + j - 1][r + k - 1].force - P.volume * P.Energy_Derivative * P.weight_gradient[i][j][k] + P.weight[i][j][k] * PID_force * P.volume;

						//Map[u + i - 1][v + j - 1][r + k - 1].force = Map[u + i - 1][v + j - 1][r + k - 1].force - P.volume * P.Energy_Derivative * P.weight_gradient[i][j][k] + P.weight[i][j][k] * (PID_force  +Dumping_force) * P.volume / cellsize[2];
						
						Max_PID_Force[0] = max(Max_PID_Force[0], PID_force[0]);
						Max_PID_Force[1] = max(Max_PID_Force[1], PID_force[1]);
						Max_PID_Force[2] = max(Max_PID_Force[2], PID_force[2]);
						if(P.force_tag = 1 && PID_force[0] != 0){
							//cout << PID_force << endl;
							Min_PID_Force[0] = min(Min_PID_Force[0], abs(PID_force[0]));
							Min_PID_Force[1] = min(Min_PID_Force[1], abs(PID_force[1]));
							Min_PID_Force[2] = min(Min_PID_Force[2], abs(PID_force[2]));
						}
						Max_Dumping_Force[0] = max(Max_Dumping_Force[0], Dumping_force[0]);
						Max_Dumping_Force[1] = max(Max_Dumping_Force[1], Dumping_force[1]);
						Max_Dumping_Force[2] = max(Max_Dumping_Force[2], Dumping_force[2]);
					}

				}
			}
		}
	}
	

}

void Grid::Compute_Velocity(){
	#pragma omp parallel for
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			for(int k = 0; k < size[2]; k++){
				if (Map[i][j][k].active && Map[i][j][k].mass > MASS_EPSILON){
					//cout << Map[i][j][k].mass << endl;
					//運動量を更新する
					Map[i][j][k].momentum = Map[i][j][k].velocity * Map[i][j][k].mass + (Map[i][j][k].force  + Map[i][j][k].g * Map[i][j][k].mass) * DT;//外力(重力項)を併せて計算
					//格子点の速度を更新する
					Map[i][j][k].velocity = Map[i][j][k].momentum / Map[i][j][k].mass;
					//cout << Map[i][j][k].velocity << endl << endl;
				}else Map[i][j][k].active = false;
		

			}
			
			

		}
	}
}

//sticky摩擦のみ適用
void Grid::Determine_Collision(){
	#pragma omp parallel for
	for (int i = 0; i < size[0]; i++){
		for (int j = 0; j < size[1]; j++){
			for(int k = 0; k < size[2]; k++){
				double y = j + Map[i][j][k].velocity[1] * DT / size[1];
				double x = i + Map[i][j][k].velocity[0] * DT / size[0];
				double z = k + Map[i][j][k].velocity[2] * DT / size[2];
				if (y < constraint || y > size[1] - constraint)
					Map[i][j][k].velocity[1] = 0;
				if (x > size[0] - constraint || x < constraint)
					Map[i][j][k].velocity[0] = 0;
				if(z < constraint || z > size[2] - constraint)
					Map[i][j][k].velocity[2] = 0;
			}

		}
	}
}


void Grid::Update_Particle_Position(){

	for (int p = 0; p < Material->points.size(); p++)
	{
		Particle& P = Material->points[p];
        double x = P.pseudo_position[0], y = P.pseudo_position[1], z = P.pseudo_position[2];
        Matrix3d affine = Matrix3d::Zero(); //アフィン速度を一時保管する変数
        int u = P.int_position[0], v = P.int_position[1], r = P.int_position[2];
        P.velocity = Vector3d::Zero();
        double density = 0;
		//#pragma omp parallel for
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
				for(int k = 0; k < 4; k++){
					if (u + i - 1 >= 0 && u + i - 1 < size[0] && v + j - 1 >= 0 && v + j - 1 < size[1] && k + r - 1 >= 0 && (k + r - 1) < size[2])  {
						if (Map[i + u - 1][j + v - 1][k + r - 1].active)
						{
							//cout << P.weight[i][j][k] << endl << endl;
							//#pragma omp critical
							P.velocity +=  Map[i + u - 1][j + v - 1][k + r - 1].velocity * P.weight[i][j][k];
							//cout << P.velocity << endl << endl;
							Eigen::Vector3d position = Eigen::Vector3d((i + u - 1)*cellsize[0] + start[0], (j + v - 1)*cellsize[1] + start[1], (k + r - 1)*cellsize[2] + start[2]);
							#pragma omp critical
							affine += P.weight[i][j][k]*Map[i + u - 1][j + v - 1][k + r - 1].velocity*(position - P.position).transpose();
						}
					}
				}

            }
        }
		P.bef_position = P.position;
		max_vel[0] = max(max_vel[0],abs(P.velocity(0)));
		max_vel[1] = max(max_vel[1],abs(P.velocity(1)));
		max_vel[2] = max(max_vel[2],abs(P.velocity(2)));
		

        P.position += P.velocity * DT;//変更

		
        P.momentum = P.velocity * P.mass;
        P.C = affine;
        P.Apply_Damping();
 


	}

}

void Grid::Update_Deformation_Gradient(){
	for (int p = 0; p < Material->points.size(); p++)
	{
		Particle& P = Material->points[p];
		double x = P.pseudo_position[0], y = P.pseudo_position[1], z = P.pseudo_position[2];
		int u = P.int_position[0], v = P.int_position[1], r = P.int_position[2];
		Matrix3d A = Matrix3d::Zero();//mapの速度勾配を足し合わせる
		#pragma omp parallel for
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				for(int k = 0; k < 4; k++){
					if (u + i - 1 >= 0 && u + i - 1 < size[0] && v + j - 1 >= 0 && v + j - 1 < size[1] && k + r - 1 >= 0 && (k + r - 1) < size[2])  {
						if (Map[i + u - 1][j + v - 1][k + r - 1].active){
							#pragma omp critical
							A += Map[i + u - 1][j + v - 1][k + r - 1].velocity * P.weight_gradient[i][j][k].transpose();//weight gradient を変える必要
						}
					}
				}

			}
		}
		P.L = A;

		A = P.L * DT + Matrix3d::Identity();
		P.F_e = A * P.F_e;
		P.J = P.F_e.determinant();


		P.volume = P.J * P.init_volume;
		P.density = P.mass / P.volume;

	}
}





/**************************************************************************************************************************/
/*****************************************************アルゴリズムの流れ******************************************************/
/**************************************************************************************************************************/

//格子点を初期化
void Grid::Init_Grid(){
    //Plastic_From_File();

    Update_From_Material();
	//getchar();
}

//格子点の状態を更新
void Grid::Update_Grid()
{
    step_g++;

	//cout << step_g << ","; //元々出していた
	logout << step_g << ",";
	//cout << sound[0] << " " << max_vel[0] << " " << max_vel[1] << " " << max_vel[2] << endl; //元々出していた
	//cout << "max pid force :"; //元々出していた
	for (int i = 0; i < 3; i++){
		// cout << " "  << std::setprecision(10) << Max_PID_Force[i] ; //元々出していた
	}
	//cout << endl; //元々出していた
	//cout << "min pid force :"; //元々出していた
	for (int i = 0; i < 3; i++){
		//cout << " "  << std::setprecision(10) << Min_PID_Force[i] ; //元々出していた
	}
	//cout << endl;
	//cout << "max dumping force :"; //元々出していた
	for (int i = 0; i < 3; i++){
		//cout << " "  << std::setprecision(10) << Max_Dumping_Force[i] ; //元々出していた
	}
	//cout << endl; //元々出していた
	if (sound[0] != max(sound[0],max(max_vel[0],max(max_vel[1],max_vel[2])))){
		cout << "before update : " <<  DT << endl;
		Calc_DT();
		//pid_dump *= 0.1;
		cout << "pid_dump : " << pid_dump << endl;
		cout << "after_update : " << DT << endl;
	}
	Calc_Error();

	//Calc_DT();
	for(int j = 0; j < 3; j++){
		max_vel[j] = 0;
		Max_PID_Force[j] = 0;
		Min_PID_Force[j] = 1e20;
		Max_Dumping_Force[j] = 0;
	}
	Update_From_Material();
	Compute_Weight_Gradient();
	Update_Force();
	Compute_Velocity();
	Determine_Collision();
	
}

//粒子の状態を更新
void Grid::Update_Particle_Cloud()
{
	
	Update_Particle_Position();
	Update_Deformation_Gradient();
	Material->Compute_Stress();
	
}


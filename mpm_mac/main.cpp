#include <iostream>
#include <cmath>
#include <cstdlib>
#include "Header/Parameters.h"
#include "Header/Grid.h"
#include <fstream>
#include <string>
#include <filesystem>
#include <sstream>
#include <omp.h>//</opt/homebrew/opt/libomp/include/omp.h>
#include <opencv2/opencv.hpp>
#include <time.h>
#include <boost/detail/workaround.hpp>//</opt/homebrew/Cellar/boost/1.80.0/include/boost/detail/workaround.hpp>
#include <boost/format.hpp>//</opt/homebrew/Cellar/boost/1.80.0/include/boost/format.hpp>
#include <numeric>
#include <filesystem>
#include <cstdio> //for error plot
#include <unistd.h>
#include <vector>
#include <algorithm>
#include <random>

#define PI 3.141592
#define OPENCV_VERSION(a,b,c) (((a) << 16) + ((b) << 8) + (c))

using namespace std;
namespace fs = std::filesystem;

bool MAX_TIFF_IMG_FLG = true;
bool OUTPUT_DATA_FLG = false;
//Particle_Cloud* object;
Particle_Cloud *object;

Grid* grid;
Particle tmp_object;
vector<Vector2d> tmp_fix;

// imageは一つ
vector<cv::Mat1w> img_target; 
// ラベルは入力を受け付ける最大臓器数による
vector<vector<cv::Mat1w> > label_target(ORGAN_KIND,vector<cv::Mat1w>()); 
//elastixにより変形させた、各labelによる画像。allも含むため臓器数＋１
vector<vector<cv::Mat1w> > elastix_img(ORGAN_KIND + 1,vector<cv::Mat1w>());
//答えに存在する、各臓器単体の画像。allも含むため臓器数＋１
vector<vector<cv::Mat1w> > answer_img(ORGAN_KIND + 1,vector<cv::Mat1w>());

//debug用
std::vector<Particle> all_points;
std::random_device rd;
std::mt19937 g(rd());       

string ii_path;
string ii_name;
string il_path;
string il_name;
string out_img;
string out_scl;
string out_csv;
string trans3D_fp;
string out3D_csv;
string outElx_img;
string ai_path;
string dice_file;
string ai_name;
string info_path;
string unconverged_path;
string divergence_path;
string init_csv;
string df_csv;
string cs_list;
string df_ratio_list;
string affine3D;
string deformationField_3D;

string target_string;

std::vector<int> time_steps;
std::vector<double> errors;
vector<double> av_time;


//画像のサイズ及びpixel spacing.Read_image_sizeが呼ばれる際に変数に値が格納される
int image_size[3] = {0,0,0};
double pixel_spacing[3] = {0,0,0};



int N = 0;
bool organ_flags[ORGAN_KIND] = {false};

double camera_x = 3;
double camera_y = 3;
double camera_z = 3;
int axis = 0;

clock_t all_time_start;
clock_t all_time_end;
clock_t one_loop_start;
clock_t one_loop_end;

//エラーのリアルタイム描画
void plotError(std::vector<int>& time_steps, std::vector<double>& errors) {
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
    if (!gnuplotPipe) {
        std::cerr << "Could not open gnuplot pipe." << std::endl;
        return;
    }
    double max_error = *std::max_element(errors.begin(), errors.end()) + 0.002;

    fprintf(gnuplotPipe, "set terminal dumb\n");  // フロントエンドをdumbに設定
    fprintf(gnuplotPipe, "set title 'Error - Time Step'\n");
    fprintf(gnuplotPipe, "set xlabel 'Time Step'\n");
    fprintf(gnuplotPipe, "set ylabel 'Error'\n");
    fprintf(gnuplotPipe, "set grid\n");
    fprintf(gnuplotPipe, "set xrange [0:%d]\n", static_cast<int>(time_steps.back())); 
    fprintf(gnuplotPipe, "set yrange [0:%.5f]\n", max_error);
    fprintf(gnuplotPipe, "plot '-' with lines title 'Error'\n");

    for (size_t i = 0; i < time_steps.size(); ++i) {
        fprintf(gnuplotPipe, "%d %f\n", time_steps[i], errors[i]);
    }

    fprintf(gnuplotPipe, "e\n");
    fflush(gnuplotPipe);
    pclose(gnuplotPipe);
}

/*
void plotError(std::vector<int>& time_steps, std::vector<double>& errors) {
    FILE* gnuplot = popen("gnuplot", "w"); 
    if (!gnuplot) {
        std::cerr << "Could not open gnuplot pipe." << std::endl;
        return;
    }

    fprintf(gnuplot, "set xrange [0:15]\n");
    fprintf(gnuplot, "set yrange [0:150]\n");
    //fprintf(gnuplot, "set terminal qt\n");

    int step = 40000;

    for (int j = 1; j < step; j += 100) {  // j=0ではlineを作れないのでj=1からとする
        fprintf(gnuplot, "plot '-' with lines\n");
        for (int i = 0; i < j && i < time_steps.size(); i += 100) {
            fprintf(gnuplot, "%d %lf\n", time_steps[i], errors[i]);
        }
        fprintf(gnuplot, "e\n"); 
        fflush(gnuplot);
        usleep(20000);
    }

    pclose(gnuplot);
}
*/
//枠描画用
int edge[][2] = {
    {0,1},
    {0,2},
    {0,3},
    {1,4},
    {1,5},
    {2,5},
    {2,7},
    {3,4},
    {3,7},
    {4,6},
    {5,6},
    {6,7}
};

GLdouble vertex[][3] = {
    {0,0,0}, //0
    {0,0,1}, //1
    {0,1,0}, //2
    {1,0,0}, //3
    {1,0,1}, //4
    {0,1,1}, //5
    {1,1,1}, //6
    {1,1,0}  //7
};


void Import_Target_Data(){

    string img_path = ii_path;
    string label_path = il_path;
    string buf,buf2;
    int flag;
    string fname;
    int k = 0;

    string img_format = ii_name;
    string label_format = il_name;

    for (k = 1; k <= image_size[2]; k++ ) {

    
        //全体を読み取る
        fname = img_path + (boost::format(img_format) % k ).str();
        img_target.push_back(cv::imread(fname,-1));       


        //label画像

        //pancreas
        if(organ_flags[0]){
            fname = label_path + "Pancreas/" + (boost::format(label_format) % k ).str();
            //fname = label_path  + (boost::format(label_format) % k ).str();
            label_target[0].push_back(cv::imread(fname,-1));       
        }
 
        //stomach
        if(organ_flags[1]){
            fname = label_path + "Stomach/" + (boost::format(label_format) % k ).str();
            label_target[1].push_back(cv::imread(fname,-1));       
        }

        //duodenum
        if(organ_flags[2]){
            fname = label_path + "Duodenum/" + (boost::format(label_format) % k ).str();
            label_target[2].push_back(cv::imread(fname,-1));       
        }
        
        if(organ_flags[3]){
            fname = label_path + "Left_Kidney/" + (boost::format(label_format) % k ).str();
            label_target[3].push_back(cv::imread(fname,-1));       
        }
        
    }
    cout  << label_target[0].size() << " " << label_target[1].size() << " " << label_target[2].size() <<  label_target[3].size() << endl;
}



//立方体の枠を描画
void Draw_area(){
    glBegin(GL_LINES);
    glColor3f(0,0,0);
    for(int i = 0; i < 12; i++){
        glVertex3dv(vertex[edge[i][0]]);
        glVertex3dv(vertex[edge[i][1]]);
    }

    glEnd();
}

//csv出力
void Output_Particle_Data_CSV(string csv){
	stringstream ss;
    ss << csv;
	ofstream pout(ss.str());
	for(int p = 0; p < object->points.size(); p++){
        pout << p << "," << object->points[p].position(0) /img_ratio[0] << "," << object->points[p].position(1) / img_ratio[1] << "," << object->points[p].position(2) / img_ratio[2] << ","<< object->points[p].z_index << "," << object->points[p].dcmval << "," << object->points[p].dcmlabel << endl;
	}
}

//Blenderなどで動画を作る際に使用する関数
//plyファイルで出力
/*
void Output_Particle_Data_ply(int n){
    stringstream ss1;
    stringstream ss2;
    stringstream ss3;
    stringstream ss4;
    stringstream ss5;
    stringstream ss6;
    string target_base = "/Users/lelab/Downloads/exam/animation";
    string stom_name = target_base + "ply_stom/output%04d.ply";
    string panc_name = target_base + "ply_panc/output%04d.ply";
    string duo_name = target_base + "ply_duo/output%04d.ply";
    string lkid_name = target_base + "ply_lkid/output%04d.ply";
    string cross_name = target_base + "ply_cross/output%04d.ply";
    string tmpcross_name = target_base + "ply_tmpcross/output%04d.ply";
    stom_name = (boost::format(stom_name)%n).str();
    panc_name = (boost::format(panc_name)%n).str();
    duo_name = (boost::format(duo_name)%n).str();
    lkid_name = (boost::format(lkid_name)%n).str();
    cross_name = (boost::format(cross_name)%n).str();
    tmpcross_name = (boost::format(tmpcross_name)%n).str();
    ss1 << panc_name;
    ofstream p1out(ss1.str());
    ss2 << stom_name;
    ofstream p2out(ss2.str());
    ss4 << tmpcross_name;
    ofstream p4out(ss4.str());
    ss5 << duo_name;
    ofstream p5out(ss5.str());
    ss6 << lkid_name;
    ofstream p6out(ss6.str());

    //cout << name << endl;

    

    //ヘッダを書き込み
    p1out << "ply" << endl;
    p1out << "format ascii 1.0" << endl;

    p2out << "ply" << endl;
    p2out << "format ascii 1.0" << endl;

    p4out << "ply" << endl;
    p4out << "format ascii 1.0" << endl;

    p5out << "ply" << endl;
    p5out << "format ascii 1.0" << endl;

    p6out << "ply" << endl;
    p6out << "format ascii 1.0" << endl;

    string vertex1 = "element vertex %d";
    string vertex2 = "element vertex %d";
    string vertex3 = "element vertex %d";
    string vertex4 = "element vertex %d";
    string vertex5 = "element vertex %d";
    string vertex6 = "element vertex %d";

    //vertex = (boost::format(vertex)%object->points.size()).str();
    vertex1 = (boost::format(vertex1)%Pancreas_num).str();
    vertex2 = (boost::format(vertex2)%Stomach_num).str();
    vertex5 = (boost::format(vertex5)%Duodenum_num).str();
    vertex6 = (boost::format(vertex6)%Left_Kidney_num).str();

    int tmp_num_data = 0;
    for(int p = 0; p < object->points.size(); p++){
        if(object->points[p].force_tag) tmp_num_data++;

    }
    
    vertex4 = (boost::format(vertex4)%tmp_num_data).str();
    p1out << vertex1 << endl;
    p2out << vertex2 << endl;
    p4out << vertex4 << endl;
    p5out << vertex5 << endl;
    p6out << vertex6 << endl;
    
    //cout << vertex << endl;
    p1out << "property float x" << endl << "property float y" << endl << "property float z" << endl;
    p1out << "property uchar red" << endl << "property uchar green" << endl << "property uchar blue" << endl;
	p1out << "end_header" << endl;

    p2out << "property float x" << endl << "property float y" << endl << "property float z" << endl;
    p2out << "property uchar red" << endl << "property uchar green" << endl << "property uchar blue" << endl;
	p2out << "end_header" << endl;

    p4out << "property float x" << endl << "property float y" << endl << "property float z" << endl;
    p4out << "property uchar red" << endl << "property uchar green" << endl << "property uchar blue" << endl;
	p4out << "end_header" << endl;

    p5out << "property float x" << endl << "property float y" << endl << "property float z" << endl;
    p5out << "property uchar red" << endl << "property uchar green" << endl << "property uchar blue" << endl;
	p5out << "end_header" << endl;


    p6out << "property float x" << endl << "property float y" << endl << "property float z" << endl;
    p6out << "property uchar red" << endl << "property uchar green" << endl << "property uchar blue" << endl;
	p6out << "end_header" << endl;

	if(N == 0){
        ss3 << cross_name;
        ofstream p3out(ss3.str());
        //制御点の数を数える
        int num_data = 0;
        for(int p = 0; p < object->points.size(); p++){
            if(object->points[p].force_tag) num_data++;

        }
        p3out << "ply" << endl;
        cout << "ply" << endl;
        p3out << "format ascii 1.0" << endl;
        vertex3 = (boost::format(vertex3)%num_data).str();
        p3out << vertex3 << endl;
        p3out << "property float x" << endl << "property float y" << endl << "property float z" << endl;
        p3out << "property uchar red" << endl << "property uchar green" << endl << "property uchar blue" << endl;
        p3out << "end_header" << endl;

        //データの読み込み
        for(int p = 0; p < object->points.size(); p++){

            
            if(object->points[p].force_tag){
                p3out << object->points[p].target_position(0) * 1 << " " << object->points[p].target_position(1) * 1 << " " << object->points[p].target_position(2) * 1 << " 0 0 0" << endl;
            }
        }
        
    }
    cout << "start output!" << endl;
    
	//データの読み込み
	for(int p = 0; p < object->points.size(); p++){
		//pout << object->points[p].position(0) * 1 << " " << object->points[p].position(1) * 1 << " " << object->points[p].position(2) * 1;
        
        if(object->points[p].dcmlabel == 1){
            p1out << object->points[p].position(0) * 1 << " " << object->points[p].position(1) * 1 << " " << object->points[p].position(2) * 1 << " 0 0 255" << endl;
        }else if(object->points[p].dcmlabel == 2){
            p2out << object->points[p].position(0) * 1 << " " << object->points[p].position(1) * 1 << " " << object->points[p].position(2) * 1 << " 255 0 0" << endl;
        }else if(object->points[p].dcmlabel == 3){
            p5out << object->points[p].position(0) * 1 << " " << object->points[p].position(1) * 1 << " " << object->points[p].position(2) * 1 << " 255 0 255" << endl;
        }else if(object->points[p].dcmlabel == 4){
            p6out << object->points[p].position(0) * 1 << " " << object->points[p].position(1) * 1 << " " << object->points[p].position(2) * 1 << " 0 255 255" << endl;
        }
        
        if(object->points[p].force_tag){
            p4out << object->points[p].position(0) * 1 << " " << object->points[p].position(1) * 1 << " " << object->points[p].position(2) * 1 << " 0 255 0" << endl;
        }
        
	}
    // exit(0);

}
*/
void Output_Particle_Data_ply(int n) {
    stringstream ss1;
    stringstream ss2;
    stringstream ss3;
    stringstream ss4;
    stringstream ss5;
    stringstream ss6;
    string target_base = "/Users/lelab/Downloads/exam/animation/";
    string stom_name = target_base + "ply_stom/output%04d.ply";
    string panc_name = target_base + "ply_panc/output%04d.ply";
    string duo_name = target_base + "ply_duo/output%04d.ply";
    string lkid_name = target_base + "ply_lkid/output%04d.ply";
    string cross_name = target_base + "ply_cross/output%04d.ply";
    string tmpcross_name = target_base + "ply_tmpcross/output%04d.ply";
    stom_name = (boost::format(stom_name) % n).str();
    panc_name = (boost::format(panc_name) % n).str();
    duo_name = (boost::format(duo_name) % n).str();
    lkid_name = (boost::format(lkid_name) % n).str();
    cross_name = (boost::format(cross_name) % n).str();
    tmpcross_name = (boost::format(tmpcross_name) % n).str();
    ss1 << panc_name;
    ofstream p1out(ss1.str());
    ss2 << stom_name;
    ofstream p2out(ss2.str());
    ss4 << tmpcross_name;
    ofstream p4out(ss4.str());
    ss5 << duo_name;
    ofstream p5out(ss5.str());
    ss6 << lkid_name;
    ofstream p6out(ss6.str());

    // Check if the file streams are open
    if (!p1out.is_open()) {
        cout << "Failed to open file: " << ss1.str() << endl;
    }
    if (!p2out.is_open()) {
        cout << "Failed to open file: " << ss2.str() << endl;
    }
    if (!p4out.is_open()) {
        cout << "Failed to open file: " << ss4.str() << endl;
    }
    if (!p5out.is_open()) {
        cout << "Failed to open file: " << ss5.str() << endl;
    }
    if (!p6out.is_open()) {
        cout << "Failed to open file: " << ss6.str() << endl;
    }

    // ここからヘッダーの書き込みとデータの処理が続きます...

    p1out << "ply" << endl;
    p1out << "format ascii 1.0" << endl;
    p2out << "ply" << endl;
    p2out << "format ascii 1.0" << endl;
    p4out << "ply" << endl;
    p4out << "format ascii 1.0" << endl;
    p5out << "ply" << endl;
    p5out << "format ascii 1.0" << endl;
    p6out << "ply" << endl;
    p6out << "format ascii 1.0" << endl;

    string vertex1 = "element vertex %d";
    string vertex2 = "element vertex %d";
    string vertex3 = "element vertex %d";
    string vertex4 = "element vertex %d";
    string vertex5 = "element vertex %d";
    string vertex6 = "element vertex %d";

    vertex1 = (boost::format(vertex1) % Pancreas_num).str();
    vertex2 = (boost::format(vertex2) % Stomach_num).str();
    vertex5 = (boost::format(vertex5) % Duodenum_num).str();
    vertex6 = (boost::format(vertex6) % Left_Kidney_num).str();

    int tmp_num_data = 0;
    for (int p = 0; p < object->points.size(); p++) {
        if (object->points[p].force_tag) tmp_num_data++;
    }

    vertex4 = (boost::format(vertex4) % tmp_num_data).str();
    p1out << vertex1 << endl;
    p2out << vertex2 << endl;
    p4out << vertex4 << endl;
    p5out << vertex5 << endl;
    p6out << vertex6 << endl;

    p1out << "property float x" << endl << "property float y" << endl << "property float z" << endl;
    p1out << "property uchar red" << endl << "property uchar green" << endl << "property uchar blue" << endl;
    p1out << "end_header" << endl;

    p2out << "property float x" << endl << "property float y" << endl << "property float z" << endl;
    p2out << "property uchar red" << endl << "property uchar green" << endl << "property uchar blue" << endl;
    p2out << "end_header" << endl;

    p4out << "property float x" << endl << "property float y" << endl << "property float z" << endl;
    p4out << "property uchar red" << endl << "property uchar green" << endl << "property uchar blue" << endl;
    p4out << "end_header" << endl;

    p5out << "property float x" << endl << "property float y" << endl << "property float z" << endl;
    p5out << "property uchar red" << endl << "property uchar green" << endl << "property uchar blue" << endl;
    p5out << "end_header" << endl;

    p6out << "property float x" << endl << "property float y" << endl << "property float z" << endl;
    p6out << "property uchar red" << endl << "property uchar green" << endl << "property uchar blue" << endl;
    p6out << "end_header" << endl;

    if (N == 0) {
        ss3 << cross_name;
        ofstream p3out(ss3.str());
        if (!p3out.is_open()) {
            cout << "Failed to open file: " << ss3.str() << endl;
        }

        int num_data = 0;
        for (int p = 0; p < object->points.size(); p++) {
            if (object->points[p].force_tag) num_data++;
        }
        p3out << "ply" << endl;
        cout << "ply" << endl;
        p3out << "format ascii 1.0" << endl;
        vertex3 = (boost::format(vertex3) % num_data).str();
        p3out << vertex3 << endl;
        p3out << "property float x" << endl << "property float y" << endl << "property float z" << endl;
        p3out << "property uchar red" << endl << "property uchar green" << endl << "property uchar blue" << endl;
        p3out << "end_header" << endl;

        for (int p = 0; p < object->points.size(); p++) {
            if (object->points[p].force_tag) {
                p3out << object->points[p].target_position(0) * 1 << " " << object->points[p].target_position(1) * 1 << " " << object->points[p].target_position(2) * 1 << " 0 0 0" << endl;
            }
        }
    }

    cout << "start output!" << endl;

    for (int p = 0; p < object->points.size(); p++) {
        if (object->points[p].dcmlabel == 1) {
            p1out << object->points[p].position(0) * 1 << " " << object->points[p].position(1) * 1 << " " << object->points[p].position(2) * 1 << " 0 0 255" << endl;
        } else if (object->points[p].dcmlabel == 2) {
            p2out << object->points[p].position(0) * 1 << " " << object->points[p].position(1) * 1 << " " << object->points[p].position(2) * 1 << " 255 0 0" << endl;
        } else if (object->points[p].dcmlabel == 3) {
            p5out << object->points[p].position(0) * 1 << " " << object->points[p].position(1) * 1 << " " << object->points[p].position(2) * 1 << " 255 0 255" << endl;
        } else if (object->points[p].dcmlabel == 4) {
            p6out << object->points[p].position(0) * 1 << " " << object->points[p].position(1) * 1 << " " << object->points[p].position(2) * 1 << " 0 255 255" << endl;
        }

        if (object->points[p].force_tag) {
            p4out << object->points[p].position(0) * 1 << " " << object->points[p].position(1) * 1 << " " << object->points[p].position(2) * 1 << " 0 255 0" << endl;
        }
    }
    // exit(0);
}
//DICE係数を計算するために、正解データを読み込む

void Read_Answer_Image(){
    string img_path = ai_path;
    string buf,buf2;
    int flag;
    string fname;
    int k = 0;
    string img_format = ai_name;
    


    for (k = 1; k <= image_size[2]; k++ ) {

        //img画像
        
        //全体を読み取る
        fname = img_path + "All/" + (boost::format(img_format) % k ).str();
        answer_img[0].push_back(cv::imread(fname,-1));       

        //pancreas
        if(organ_flags[0]){
            fname = img_path + "Pancreas/" + (boost::format(img_format) % k ).str();
            //fname = img_path  + (boost::format(img_format) % k ).str();
            answer_img[1].push_back(cv::imread(fname,-1));       
        }
 
        //stomach
        if(organ_flags[1]){
            fname = img_path + "Stomach/" + (boost::format(img_format) % k ).str();
            answer_img[2].push_back(cv::imread(fname,-1));       
        }

        //duodenum
        if(organ_flags[2]){
            fname = img_path + "Duodenum/" + (boost::format(img_format) % k ).str();
            answer_img[3].push_back(cv::imread(fname,-1));       
        }
        //left_kidney
        if(organ_flags[3]){
            fname = img_path + "Left_Kidney/" + (boost::format(img_format) % k ).str();
            answer_img[4].push_back(cv::imread(fname,-1));       
        }
    }
    cout  << answer_img[0].size() << " " << answer_img[1].size() << " " << answer_img[2].size() << " " << answer_img[3].size() << " " << answer_img[4].size() << endl;
}

double calc_dice(vector<cv::Mat1w> img1, vector<cv::Mat1w> img2, string title){
    double same_region = 0;
    double all_region = 0;
    cout << "calc_dice " << endl;
    //cout << img1.size() << " " << img1[0].cols << " " << img1[0].rows << endl;
     for (int k = 0; k < img1.size(); k++){
        for(int i = 0; i < img1[k].cols; i++){
            for(int j = 0; j < img1[k].rows; j++){
                if (img1[k].at<unsigned short>(j,i)!= 0 && img2[k].at<unsigned short>(j,i) != 0){
                    same_region++;
                }
                if(img1[k].at<unsigned short>(j,i) != 0 || img2[k].at<unsigned short>(j,i) != 0){
                    all_region++;
                }
            }
        }
    }
    all_region += same_region;
    cout << "----------" << title <<  "----------" << endl;
    logout << "----------" << title <<  "----------" << endl;
    cout << "same_region " << same_region << " all_region " << all_region << endl;
    logout << "same_region " << same_region << " all_region " << all_region << endl;
    cout << "Dice Score is " << 2*same_region / all_region << endl;
    logout << "Dice Score is " << 2*same_region / all_region << endl;
    return 2*same_region / all_region;

}

//elastixによる3次元画像変換後の情報を出力
void Calc_3D_Trans_Image(){
	MatrixXd warp(4,4);
	MatrixXd trans(4,4);
	MatrixXd center(4,4);
	MatrixXd m_center(4,4);
	MatrixXd affine(4,4);

	MatrixXd preprocess(4,4);

    vector<double> sp(2,0);
	vector<double> cr(2,0);

	vector<double> num;
    bool tr_flag = true;

	string elastix3d_path = trans3D_fp;
	cout << "elastix3dpath : " << elastix3d_path << endl;
	string buf,buf2, param;
	ifstream read_file(elastix3d_path);
	cout << "affine3d" << endl;
	while(getline(read_file,buf)){

		if(buf[0] == '('){
			buf.erase(0,1);
			buf.erase(buf.size()-1);
		}else if(buf[0] == '/'){
			buf.erase();
		}
		num.resize(12,0);

		if(buf.find("TransformParameters") != std::string::npos && tr_flag){
			istringstream iss(buf);
			iss >> param >> num[0] >> num[1] >> num[2] >> num[3] >> num[4] >> num[5] >> num[6] >> num[7] >> num[8] >> num[9] >> num[10] >> num[11];
			if (param == "TransformParameters"){
				cout << param << " " << num[0] << " " << num[1] << " " << num[2] << " " << num[3] << " " << num[4] << " " << num[5] << " " << num[6] << " " << num[7] << " " << num[8] << " " << num[9] << " " << num[10] << " " << num[11] << endl;
				logout << param << " " << num[0] << " " << num[1] << " " << num[2] << " " << num[3] << " " << num[4] << " " << num[5] << " " << num[6] << " " << num[7] << " " << num[8] << " " << num[9] << " " << num[10] << " " << num[11] << endl;
				tr_flag = false;
			}
		}
		if(buf.find("Spacing") != std::string::npos){
			istringstream iss(buf);
			iss >> param >> sp[0] >> sp[1] >> sp[2];
			//cout << param << sp[0] << sp[1] << endl;
			if (param == "Spacing"){
				cout << param << " " << sp[0] << " " << sp[1] << " " << sp[2] << endl;
				logout << param << " " << sp[0] << " " << sp[1] << " " << sp[2]<< endl;
				//getchar();
			}
		}
		if(buf.find("CenterOfRotationPoint") != std::string::npos){
			istringstream iss(buf);
			iss >> param >> cr[0] >> cr[1] >> cr[2];
			//cout << param << cr[0] << cr[1] << endl;
			if (param == "CenterOfRotationPoint"){
				cout << param << " " << cr[0] << " " << cr[1] << " " << cr[2] << endl;
				logout << param << " " << cr[0] << " " << cr[1] << " " << cr[2]<< endl;
				//getchar();
				break;
			}
			
		}
			
	}

	affine << num[0] ,num[1], num[2], 0,
			num[3], num[4], num[5], 0,
			num[6], num[7], num[8], 0,
			0, 0, 0, 1;
	trans << 1, 0, 0, num[9],
			0, 1, 0, num[10],
			0, 0, 1, num[11],
			0, 0, 0, 1;



	VectorXd rotate_center(4);
    VectorXd spacing(3);
	VectorXd vec = VectorXd::Zero(4);
	VectorXd new_vec = VectorXd::Zero(4);
	VectorXd new_trans(4);
	
    spacing << sp[0], sp[1], sp[2];
	rotate_center << cr[0], cr[1], cr[2], 0;

    //どの場合によっても共通の行列
	preprocess << 1, 0, 0, 0,
				0, -1, 0, 0,
				0, 0, 1, 0,
				0, 0, 0, 1;
	warp = preprocess * affine * preprocess;
	trans = preprocess * trans * preprocess;
	new_trans << trans(0,3), trans(1,3), trans(2,3),1;
    //出力用のデータを作成

	stringstream ss;
	ss << out3D_csv;
    ofstream pout(ss.str());
    for(int i = 0; i < elastix_img.size(); i++){
        for(int j = 0; j < image_size[2]; j++){
            elastix_img[i].push_back(cv::Mat::zeros(cv::Size(image_size[0],image_size[1]),CV_16U));
        }
    }
	for(int p = 0; p < object->points.size(); p++){
		Particle& P = object->points[p];

		vec(0) = P.init_position[0] / img_ratio[0] * spacing(0);
		vec(1) = P.init_position[1] / img_ratio[1] * spacing(1);
		vec(2) = P.init_position[2] / img_ratio[2] * spacing(2);
		vec(3) = 1;
		vector<double> out_pos(3,0);
		new_vec = vec - rotate_center;
		new_vec = warp * new_vec;
		new_vec += rotate_center;
		
		new_vec += new_trans;
		out_pos[0] = new_vec(0) /(spacing(0));
		out_pos[1] = new_vec(1) /(spacing(1));
		out_pos[2] = new_vec(2) / spacing(2);

		//画像出力用
		int i = new_vec(0) / spacing(0);
		int j = new_vec(1) / spacing(1); 
		int z = new_vec(2) / spacing(2);
		unsigned short val = P.dcmval;   
		elastix_img[0][z].at<unsigned short>(image_size[1]-j,i) = 65535; 
        elastix_img[P.dcmlabel][z].at<unsigned short>(image_size[1]-j,i) = 65535; 
		pout << p << "," << out_pos[0] << "," << out_pos[1] << "," << out_pos[2] << "," << P.z_index << "," << val << "," << P.dcmlabel << endl;    
	}

    //ファイル書き込み用
    ofstream writing_file;
    //ファイルが存在するか確認し、存在しなければ作成
    cout << "dice" << fs::exists(dice_file) << endl;
    if(!fs::exists(dice_file)){
        writing_file.open(dice_file,std::ios::out);
    }else{
        writing_file.open(dice_file,std::ios::app);
    }

    vector<string> elx_titles{"elastix-ans:all","elastix-ans:pancreas","elastix-ans:stomach","elastix-ans:duodenum","elastix-ans:left_kidney"};
    vector<double> dice_score(ORGAN_KIND + 1,0);

    for(int i = 0; i < dice_score.size(); i++){
        dice_score[i] = calc_dice(elastix_img[i],answer_img[i], elx_titles[i]);
        cout << dice_score[i] << endl;
    }
    string title_name; 

    for (int j = 0; j < dice_score.size(); j++){
        title_name = elx_titles[j];
        writing_file << to_string(cs) << "," << target_string << "," << title_name << "," << to_string(dice_score[j]) << endl;
        cout << to_string(cs) << "," << target_string << "," << title_name << "," << to_string(dice_score[j]) << endl;
    }
    //writing_file << endl;

    writing_file.close();

    //　ダイス係数を計算できるように、臓器ごとのデータは本来の画素値ではないが、outimg先に保存する
    int err = 0;
    string dirname = outElx_img;
    string fpath;
    string filename_base = "%04d.tiff";
    string filename;


    system(("mkdir " + dirname + "All").c_str());
    //ファイル名を作成
    fpath = dirname + "/All/";

    for(int l = 1; l <= image_size[2]; l++){
        filename = fpath  + (boost::format(filename_base) % l ).str();
        cv::imwrite(filename,elastix_img[0][l-1]);
    }

    // 臓器ごとのフォルダを作成＆データ保存
    if(organ_flags[0]){

        system(("mkdir " + dirname + "Pancreas").c_str());
        //ファイル名を作成
        fpath = dirname + "/Pancreas/";
        for(int l = 1; l <= image_size[2]; l++){
            filename = fpath  + (boost::format(filename_base) % l ).str();
            cv::imwrite(filename,elastix_img[1][l-1]);
        }
    }
    if(organ_flags[1]){
        system(("mkdir " + dirname + "Stomach").c_str());
        //ファイル名を作成
        fpath = dirname + "/Stomach/";
        for(int l = 1; l <= image_size[2]; l++){
            filename = fpath  + (boost::format(filename_base) % l ).str();
            cv::imwrite(filename,elastix_img[2][l-1]);
        }
    }
    if(organ_flags[2]){
        system(("mkdir " + dirname + "Duodenum").c_str());
        //ファイル名を作成
        fpath = dirname + "/Duodenum/";
        for(int l = 1; l <= image_size[2]; l++){
            filename = fpath  + (boost::format(filename_base) % l ).str();
            cv::imwrite(filename,elastix_img[3][l-1]);
        }
    }
    if(organ_flags[3]){
        system(("mkdir " + dirname + "Left_Kidney").c_str());
        //ファイル名を作成
        fpath = dirname + "/Left_Kidney/";
        for(int l = 1; l <= image_size[2]; l++){
            filename = fpath  + (boost::format(filename_base) % l ).str();
            cv::imwrite(filename,elastix_img[4][l-1]);
        }
    }
    cout << elastix_img.size() << " " << elastix_img[0].size() << endl;

	cout << "calc_end" << endl;
    return;

}


void Output_Particle_Data(){

    vector<vector<cv::Mat1w> > img(ORGAN_KIND + 1,vector<cv::Mat1w>());
    vector<vector<cv::Mat1w> > rescale(ORGAN_KIND + 1,vector<cv::Mat1w>());



    for(int j = 0; j < ORGAN_KIND + 1; j++){
        for(int i = 0; i < image_size[2]; i++){
            img[j].push_back(cv::Mat::zeros(cv::Size(image_size[0],image_size[1]),CV_16U));
            rescale[j].push_back(cv::Mat::zeros(cv::Size(image_size[0],image_size[1]),CV_16U));
        }
    }
    
    int particle_num = 0;
    for(int l = 0; l < object->points.size(); l++){

        int i = object->points[l].position[0] / img_ratio[0];
        int j = object->points[l].position[1] / img_ratio[1];
        int z_index = object->points[l].position[2] / img_ratio[2];
        unsigned short val = object->points[l].dcmval; 
        int id = object->points[l].dcmlabel;
        //cout << i << " " << j << " "<< z_index << endl;
        img[0][z_index].at<unsigned short>(img[0][z_index].rows-j,i) = val; 
        rescale[0][z_index].at<unsigned short>(rescale[0][z_index].rows-j,i) = 65535;
        if(val > 0)particle_num++;    
        img[id][z_index].at<unsigned short>(img[0][z_index].rows-j,i) = val;
        rescale[id][z_index].at<unsigned short>(rescale[0][z_index].rows-j,i) = 65535;

    }

    //dice係数を計算する
    //セグメンテーションで画素値がない部分にもラベルがある場合が存在するため、
    //rescaleを使用している.そのため読み込みもラベル、
    //もしくはそれに準ずる画像である必要がある
    if(mode == 2){
        //ファイル書き込み用
        ofstream writing_file;
        //ファイルが存在するか確認し、存在しなければ作成
        cout << "dice" << fs::exists(dice_file) << endl;
        if(!fs::exists(dice_file)){
            writing_file.open(dice_file,std::ios::out);
        }else{
            writing_file.open(dice_file,std::ios::app);
        }

        vector<string> sim_titles{"sim-ans:all","sim-ans:pancreas","sim-ans:stomach","sim-ans:duodenum","sim-ans:left_kidney"};
        vector<double> dice_score(ORGAN_KIND + 1,0);

        dice_score[0] = calc_dice(rescale[0],answer_img[0], sim_titles[0]);
        if(organ_flags[0])dice_score[1] = calc_dice(rescale[1],answer_img[1], sim_titles[1]);
        if(organ_flags[1])dice_score[2] = calc_dice(rescale[2],answer_img[2], sim_titles[2]);
        if(organ_flags[2])dice_score[3] = calc_dice(rescale[3],answer_img[3], sim_titles[3]);
        if(organ_flags[3])dice_score[4] = calc_dice(rescale[4],answer_img[4], sim_titles[4]);

        string title_name; 
        for (int j = 0; j < dice_score.size(); j++){
            title_name = sim_titles[j];
            writing_file << to_string(cs) << "," << target_string << "," << title_name << "," << to_string(dice_score[j]) << endl;
        }


        writing_file.close();
    }
    else if(mode == 4){
        //ファイル書き込み用
        ofstream writing_file;
        //ファイルが存在するか確認し、存在しなければ作成
        cout << "dice" << fs::exists(dice_file) << endl;
        if(!fs::exists(dice_file)){
            writing_file.open(dice_file,std::ios::out);
        }else{
            writing_file.open(dice_file,std::ios::app);
        }

        vector<string> sim_titles{"sim-ans:all","sim-ans:pancreas","sim-ans:stomach","sim-ans:duodenum","sim-ans:left_kidney"};
        vector<double> dice_score(ORGAN_KIND + 1,0);

        dice_score[0] = calc_dice(rescale[0],answer_img[0], sim_titles[0]);
        if(organ_flags[0])dice_score[1] = calc_dice(rescale[1],answer_img[1], sim_titles[1]);
        if(organ_flags[1])dice_score[2] = calc_dice(rescale[2],answer_img[2], sim_titles[2]);
        if(organ_flags[2])dice_score[3] = calc_dice(rescale[3],answer_img[3], sim_titles[3]);
        if(organ_flags[3])dice_score[4] = calc_dice(rescale[4],answer_img[4], sim_titles[4]);

        string title_name; 
        /*
        for (int j = 0; j < dice_score.size(); j++){
            title_name = sim_titles[j];
            writing_file << to_string(cs) << "," << target_string << "," << title_name << "," << to_string(dice_score[j]) << endl;
        }
        */
        writing_file << to_string(cs_x) << "," << to_string(cs_y) << "," << to_string(cs_z) << "," << to_string(dice_score[1]) << endl;



        writing_file.close();
    }

    if (OUTPUT_DATA_FLG){
        //画像出力ステップ
        string fpath;
        string filename_base = "%04d.tiff";
        string filename;
        string dirname;
        int err = 0;

        for (int i = 0; i < 2; i++){
            //MAX_TIFFフラグが立っていない時はスキップ
            if(i == 0) dirname = out_img;
            else if(!MAX_TIFF_IMG_FLG && i == 1) continue; 
            else dirname = out_scl;
            //最初のステップだった場合、フォルダ内を初期化する
            err = 0;
            if(mode == 1)dirname = dirname + "Answer";
            if(mode == 2)dirname = dirname + "Estimation";
    
            string command1 = ("rm -rf " + dirname);
            err = system(command1.c_str());
            if(err != 0){
                cerr << "ディレクトリを削除できませんでした" << endl;
                exit(1);
            }

            string command2 = ("mkdir " + dirname);
            err = system(command2.c_str());
            if(err != 0){
                cerr << "ディレクトリを作成出来ませんでした" << endl;
                exit(1);
            }
            
            
            system(("mkdir " + dirname + "/All").c_str());
            //ファイル名を作成
            string fpath = dirname + "/All/";

            for(int l = 1; l <= image_size[2]; l++){
                filename = fpath  + (boost::format(filename_base) % l ).str();
                if(i == 0) cv::imwrite(filename,img[0][l-1]);
                else if(i == 1) cv::imwrite(filename,rescale[0][l-1]);
            }

            // 臓器ごとのフォルダを作成＆データ保存
            if(organ_flags[0]){

                system(("mkdir " + dirname + "/Pancreas").c_str());
                //ファイル名を作成
                fpath = dirname + "/Pancreas/";
                for(int l = 1; l <= image_size[2]; l++){
                    filename = fpath  + (boost::format(filename_base) % l ).str();
                    if(i == 0) cv::imwrite(filename,img[1][l-1]);
                    else if(i == 1) cv::imwrite(filename,rescale[1][l-1]);
                }
            }
            if(organ_flags[1]){
                system(("mkdir " + dirname + "/Stomach").c_str());
                //ファイル名を作成
                fpath = dirname + "/Stomach/";
                for(int l = 1; l <= image_size[2]; l++){
                    filename = fpath  + (boost::format(filename_base) % l ).str();
                    if(i == 0) cv::imwrite(filename,img[2][l-1]);
                    else if(i == 1) cv::imwrite(filename,rescale[2][l-1]);
                }
            }
            if(organ_flags[2]){
                system(("mkdir " + dirname + "/Duodenum").c_str());
                //ファイル名を作成
                fpath = dirname + "/Duodenum/";
                for(int l = 1; l <= image_size[2]; l++){
                    filename = fpath  + (boost::format(filename_base) % l ).str();
                    if(i == 0) cv::imwrite(filename,img[3][l-1]);
                    else if(i == 1) cv::imwrite(filename,rescale[3][l-1]);
                }
            }
            if(organ_flags[3]){
                system(("mkdir " + dirname + "/Left_Kidney").c_str());
                //ファイル名を作成
                fpath = dirname + "/Left_Kidney/";
                for(int l = 1; l <= image_size[2]; l++){
                    filename = fpath  + (boost::format(filename_base) % l ).str();
                    if(i == 0) cv::imwrite(filename,img[4][l-1]);
                    else if(i == 1) cv::imwrite(filename,rescale[4][l-1]);
                }
            }
        }
    }   
    Output_Particle_Data_CSV(out_csv);
}


void keybord(unsigned char key, int x, int y){
    switch(key) {
        case 'q':
        case 'Q':
        case '\033':
            exit(0);
            break;
        case 'x':
            axis = 0;
            break;
        case 'y':
            axis = 1;
            break;
        case 'z':
            axis = 2;
            break;
        case 'p':
            getchar();
            break;
        case 'o':
            cout << "output!" << endl;
            //Output_Particle_Data_CSV();
            break;
        default :
            break;
    }
}

void mouse(int button, int state, int x, int y){
    if((button == GLUT_LEFT_BUTTON) && (state == GLUT_DOWN)){
        if(axis == 0) camera_x += 0.5;
        if(axis == 1) camera_y += 0.5;
        if(axis == 2) camera_z += 0.5;
    }

    if((button == GLUT_RIGHT_BUTTON) && (state == GLUT_DOWN)){
        if(axis == 0) camera_x -= 0.5;
        if(axis == 1) camera_y -= 0.5;
        if(axis == 2) camera_z -= 0.5;
    }
  
}

void Calc_Init_Error(){
    //ファイル書き込み用
    ofstream writing_file;
    //ファイルが存在するか確認し、存在しなければ作成
    cout << "dice" << fs::exists(dice_file) << endl;
    if(!fs::exists(dice_file)){
        writing_file.open(dice_file,std::ios::out);
    }else{
        writing_file.open(dice_file,std::ios::app);
    }

    vector<string> init_titles{"init-ans:pancreas","init-ans:stomach","init-ans:duodenum","init-ans:left_kidney"};
    vector<double> dice_score(ORGAN_KIND,0);

    for(int i = 0; i < dice_score.size(); i++){
        dice_score[i] = calc_dice(label_target[i],answer_img[i+1], init_titles[i]);
        cout << dice_score[i] << endl;
    }
    string title_name; 

    for (int j = 0; j < dice_score.size(); j++){
        title_name = init_titles[j];
        writing_file << to_string(cs) << "," << target_string << "," << title_name << "," << to_string(dice_score[j]) << endl;
        cout << to_string(cs) << "," << target_string << "," << title_name << "," << to_string(dice_score[j]) << endl;
    }
    //writing_file << endl;

    writing_file.close();

    

}

//https://suzulang.com/save-to-ppm-by-glreadpixels/

//画像の上下を反転する
void reverse_Y(const int width, const int height, unsigned char* p) {

  int WidthByte = width * 3;

  size_t HalfHEIGHT = height / 2;
  for (int ha = 0; ha < HalfHEIGHT; ha++) {

    int hb = (height - ha) - 1;

    unsigned char* pha = p + ha * WidthByte;
    unsigned char* phb = p + hb * WidthByte;

    if (ha != hb) {
      for (size_t i = 0; i < WidthByte; i++) {
        std::swap(pha[i], phb[i]);
      }
    }

  }

}
//画像描画用
void pnmP3_Write(const char* const fname, const int vmax, const int width, const int height, const unsigned char* const p) { // PPM ASCII


  FILE* fp = fopen(fname, "w");

  fprintf(fp, "P3\n%d %d\n%d\n", width, height, vmax);

  size_t k = 0;
  for (size_t i = 0; i < (size_t)height; i++) {
    for (size_t j = 0; j < (size_t)width; j++) {
      fprintf(fp, "%d %d %d ", p[k * 3 + 0], p[k * 3 + 1], p[k * 3 + 2]);
      k++;
    }
    fprintf(fp, "\n");
  }

  fclose(fp);


}

void saveGLtoPPM(const int width, const int height, const char* pathname) {

  //メモリ確保
  std::vector<GLubyte> buffer;
  buffer.resize(width*height * 3);

  //画像取得
  glReadBuffer(GL_FRONT);
  glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, buffer.data());

  //上下反転
  reverse_Y(width, height, buffer.data());

  pnmP3_Write(
    pathname,
    255,
    width,
    height,
    buffer.data()); // PPM ASCII

}

//GLの初期設定
void InitGL(){
    int dummy = 0;
    glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
    glutInitWindowPosition(0,0);
    glutInit(&dummy,NULL);
    glutCreateWindow("cross-sectional driven simulation");
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE | GLUT_DEPTH);
    glClearColor(1,1,1,0);
    glEnable(GL_DEPTH_TEST);
    glutKeyboardFunc(keybord);
    glutMouseFunc(mouse);
    glEnable(GL_POINT_SMOOTH);
    glPointSize(3);
    glMatrixMode(GL_PROJECTION);
    cout << "InitGL complete!" << endl;
}



//レンダリング関数
void Render(){
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glViewport(0,0,WINDOW_WIDTH, WINDOW_HEIGHT);
    glLoadIdentity();
    //glOrtho(0, VIEW_WIDTH, 0, VIEW_HEIGHT, 0, 1); //最初
    //glOrtho(-2.0, 2.0, -2.0, 2.0, -2.0, 2.0);
    gluPerspective(30.0, (double)VIEW_WIDTH / (double)VIEW_HEIGHT, 1.0, 100.0);
    //gluLookAt(camera_x + 2, camera_y-5.5, camera_z - 2.5, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    gluLookAt(camera_x + 2, -(camera_y-8.5), camera_z - 2.5, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    //gluLookAt(camera_x, camera_y, camera_z, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    //gluLookAt(0.5, 0.5, 5, 0.5, 0.5, 0.5, 0.0, 1.0, 0.0);
    //gluLookAt(5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 1.0);
    
    //物体を描画
    //cout << "object: " << object << endl;
    object->Draw();
    Draw_area();
    //拘束、境界を描画
    //fix->Draw();
    //getchar();
    
    glFlush();

    GLenum err;
    while ((err = glGetError()) != GL_NO_ERROR) {
        std::cerr << "OpenGLエラー: " << err << std::endl;
    }
}


//プログラムを終了させる関数
void terminate_program(){
    cout << "time step : " << DT << endl;
    logout << "time step : " << DT << endl;

    double av = std::accumulate(av_time.begin(),av_time.end(),0.0) / av_time.size();
    cout << "average time : " << av << "sec." << endl;
    logout << "average time : " << av << "sec." << endl;
    cout << "total time : " << (double)(all_time_end - all_time_start) / CLOCKS_PER_SEC << "sec." << endl;
    logout << "total time : " << (double)(all_time_end - all_time_start) / CLOCKS_PER_SEC << "sec." << endl;
    exit(0);
}

void Update(){
    one_loop_start = clock();
    //Output_Particle_Data_ply(N);
    /*
    Output_Particle_Data_CSV();
    cout << "out_csv" << endl;
    exit(0);
    */
    if(mode == 2 && N == 0){
        saveGLtoPPM(WINDOW_WIDTH, WINDOW_HEIGHT, "start.ppm");
    }
    if (mode == 1 && N == 0){
        cout << "calc init position" << endl;
        Output_Particle_Data_CSV(init_csv);
    }
    if (mode == 3){
        cout << "calc init error" << endl;
        Calc_Init_Error();
        cout << "calc_elastix_data" << endl;
        Calc_3D_Trans_Image();
        exit(0);
    }
    if(pid_output_flag){
        //saveGLtoPPM(WINDOW_WIDTH, WINDOW_HEIGHT, "end.ppm");
        cout << "output! N= " << N << endl;
        // elastixによる変形データを作成する
        if(mode==2 || mode == 4){
            //unconverged_flag = true;
            if(unconverged_flag){
                logout << "write unconverged info" << endl;
                cout << "write unconverged info" << endl;
                //ファイル書き込み用
                ofstream unconverged_file;
                //ファイルが存在するか確認し、存在しなければ作成
                cout << "dice" << fs::exists(unconverged_path) << endl;
                if(!fs::exists(unconverged_path)){
                    unconverged_file.open(unconverged_path,std::ios::out);
                }else{
                    unconverged_file.open(unconverged_path,std::ios::app);
                }
                unconverged_file << to_string(cs);
                if(organ_flags[0]) unconverged_file << ",Pancreas";
                if(organ_flags[1]) unconverged_file << ",Stomach";
                if(organ_flags[2]) unconverged_file << ",Duodenum";
                if(organ_flags[3]) unconverged_file << ",Left_Kidney";
                unconverged_file << endl;
                unconverged_file.close();
            }
            if(divergence_flag){
                logout << "write divergence info" << endl;
                cout << "write divergence info" << endl;
                //ファイル書き込み用
                ofstream divergence_file;
                //ファイルが存在するか確認し、存在しなければ作成
                cout << "dice" << fs::exists(unconverged_path) << endl;
                if(!fs::exists(unconverged_path)){
                    divergence_file.open(divergence_path,std::ios::out);
                }else{
                    divergence_file.open(divergence_path,std::ios::app);
                }
                divergence_file << to_string(cs);
                if(organ_flags[0]) divergence_file << ",Pancreas";
                if(organ_flags[1]) divergence_file << ",Stomach";
                if(organ_flags[2]) divergence_file << ",Duodenum";
                if(organ_flags[3]) divergence_file << ",Left_Kidney";
                divergence_file << endl;
                divergence_file.close();
            }
            cout << "Import Answer Data !" << endl;
            Read_Answer_Image();
            //cout << "Calculate elastix output !" << endl;
            //Calc_3D_Trans_Image();
        }

        Output_Particle_Data();
        

        terminate_program();

    }

    //MPMアルゴリズムに沿って更新
    grid->Update_Grid();
    grid->Update_Particle_Cloud();
    N++;
    one_loop_end = clock();
    all_time_end = clock();
    av_time.push_back((double)(one_loop_end - one_loop_start) / CLOCKS_PER_SEC);
    time_steps.push_back(N);
    errors.push_back(grid->error);
    plotError(time_steps, errors);

    cout << "N= " << N << endl;
    //glutPostRedisplay();
}

//離散化の際に使う倍率を決定する
void CalcImgRatio(){
    double div = max(max(image_size[0]*pixel_spacing[0],image_size[1]*pixel_spacing[1]),image_size[2]*pixel_spacing[2]);
    img_ratio[0] = pixel_spacing[0] / div;
    img_ratio[1] = pixel_spacing[1] / div;
    img_ratio[2] = pixel_spacing[2] / div;
    cout << "img_ratio[0], img_ratio[1], img_ratio[2]" << img_ratio[0] << ", " << img_ratio[1] << ", " << img_ratio[2] << endl;
}
//MPMの情報を初期化。ここで物体の離散化を行う
void InitMPM(){
    object = new Particle_Cloud();

    //////////////////////////////////////////////////////////
    /*                      弾性体を離散化                      */
    //////////////////////////////////////////////////////////
   
    Import_Target_Data();
    //画像の倍率を計算
    CalcImgRatio();
    
    Vector3d vel;
    
    vel << 0,0,0;

    

    int import_label = 0;

    //Pancreas = 1
    //Stomach = 2
    //Duodenum = 3
    // 0.7895 : 2.4 = 1/576 : 1 / x
    import_label = 0;
    int label = 0;
    enum MaterialType organ_type;
    int min_i = 10000;
    int min_j = 10000;
    int max_i = 0;
    int max_j = 0;

    for(int n = 0; n < ORGAN_KIND; n++){
        if(!organ_flags[n])continue;
        if(n == 0) organ_type = PANC;
        else if(n == 1) organ_type = STOM;
        else if(n == 2) organ_type = DUO;
        else if(n == 3) organ_type = LKID;
        //for(int k = 0; k < image_size[2]; k++){
        for(int k = 0; k < image_size[2]; k++){
            for(int i = 0; i < label_target[n][k].cols; i++){
                for(int j = 0; j < label_target[n][k].rows; j++){
                    if(label_target[n][k].at<unsigned short>(j,i) != 0){
                        
                        if(import_label == 0){
                            if(n == 0){
                                min_i = min(min_i, i);
                                min_j = min(min_j, j);
                                max_i = max(max_i, i);
                                max_j = max(max_j, j);
                            }
                                
                            if (i == cs_x || label_target[n][k].rows - j == cs_y || k == cs_z){
                                tmp_object = Particle(Vector3d(img_ratio[0] * i, img_ratio[1] * (label_target[n][k].rows - j),img_ratio[2] * k), vel, 1, 1, (n+1), organ_type,img_target[k].at<unsigned short>(j,i),k+1, true);
                            }else{
                                tmp_object = Particle(Vector3d(img_ratio[0] * i, img_ratio[1] * (label_target[n][k].rows - j),img_ratio[2] * k), vel, 1, 1, (n+1), organ_type,img_target[k].at<unsigned short>(j,i),k+1, true);
                            }
                            
                            //tmp_object = Particle(Vector3d(i, (label_target[n][k].rows - j), k), vel, 1, 1, (n+1), organ_type,img_target[k].at<unsigned short>(j,i),k+1);
                            object->points.push_back(tmp_object);
                            all_points.push_back(tmp_object);
                            if(n == 0) Pancreas_num++;
                            else if (n == 1) Stomach_num++;
                            else if(n == 2) Duodenum_num++;
                            else if(n == 3) Left_Kidney_num++;
                            //import_label++;
                        }else if (import_label < 10){
                            import_label++;
                        }else{
                            import_label = 0;
                        }

                    }
                }
            }
        }
    }

    /*デバッグ用！！objectから  ランダムにn個のだけのオブジェクトを作る
    1264s行目、フルでオブジェクトを作る時は   object->points.push_back(tmp_object);
    ランダムな数のオブジェクトを作るときは      all_points.push_back(tmp_object);
    */
    std::shuffle(all_points.begin(), all_points.end(), g);
    int num_points_to_add = std::min(1000000, (int)all_points.size());

    for(int i = 0; i < num_points_to_add; i++) {
        object->points.push_back(all_points[i]);
    }
    //randomオブジェクト生成ここまで

    cout << "min_i : " << min_i << ", max_i " << max_i << endl;
    cout << "min_j : " << min_j << ", max_j " << max_j << endl;

    cout << "Stomach_num " << Stomach_num << endl;
    cout << "Pancreas_num " << Pancreas_num << endl;
    cout << "Duodenum_num " << Duodenum_num << endl;
    cout << "Left_Kidney_num " << Left_Kidney_num << endl;
    logout << "Stomach_num " << Stomach_num << endl;
    logout << "Pancreas_num " << Pancreas_num << endl;
    logout << "Duodenum_num " << Duodenum_num << endl;
    logout << "Left_Kidney_num " << Left_Kidney_num << endl;
    //getchar();



    //格子点を初期化
    grid = new Grid(Vector3d(0,0,0), Vector3d(VIEW_WIDTH, VIEW_HEIGHT, 1), Vector3i(121,121,121), object);
    grid->Init_Grid();


    vector<vector<cv::Mat1w> > img(ORGAN_KIND + 1,vector<cv::Mat1w>());
    vector<vector<cv::Mat1w> > rescale(ORGAN_KIND + 1,vector<cv::Mat1w>());


    Read_Answer_Image();

    for(int j = 0; j < ORGAN_KIND + 1; j++){
        for(int i = 0; i < image_size[2]; i++){
            img[j].push_back(cv::Mat::zeros(cv::Size(image_size[0],image_size[1]),CV_16U));
            rescale[j].push_back(cv::Mat::zeros(cv::Size(image_size[0],image_size[1]),CV_16U));
        }
    }
    
    int particle_num = 0;
    for(int l = 0; l < object->points.size(); l++){

        int i = object->points[l].target_position[0] / img_ratio[0];
        int j = object->points[l].target_position[1] / img_ratio[1];
        int z_index = object->points[l].target_position[2] / img_ratio[2];
        unsigned short val = object->points[l].dcmval; 
        int id = object->points[l].dcmlabel;
        //cout << i << " " << j << " "<< z_index << endl;
        img[0][z_index].at<unsigned short>(img[0][z_index].rows-j,i) = val; 
        rescale[0][z_index].at<unsigned short>(rescale[0][z_index].rows-j,i) = 65535;
        if(val > 0)particle_num++;    
        img[id][z_index].at<unsigned short>(img[0][z_index].rows-j,i) = val;
        rescale[id][z_index].at<unsigned short>(rescale[0][z_index].rows-j,i) = 65535;

    }

    vector<string> tar_titles{"tar-ans:all","tar-ans:pancreas","tar-ans:stomach","tar-ans:duodenum","tar-ans:left_kidney"};
    vector<double> dice_score(ORGAN_KIND + 1,0);

    dice_score[0] = calc_dice(rescale[0],answer_img[0], tar_titles[0]);
    if(organ_flags[0])dice_score[1] = calc_dice(rescale[1],answer_img[1], tar_titles[1]);
    if(organ_flags[1])dice_score[2] = calc_dice(rescale[2],answer_img[2], tar_titles[2]);
    if(organ_flags[2])dice_score[3] = calc_dice(rescale[3],answer_img[3], tar_titles[3]);
    if(organ_flags[3])dice_score[4] = calc_dice(rescale[4],answer_img[4], tar_titles[4]);

    for (int j = 0; j < dice_score.size(); j++){
        string title_name;
        title_name = tar_titles[j];
        //writing_file << to_string(cs) << "," << target_string << "," << title_name << "," << to_string(dice_score[j]) << endl;
        cout << to_string(cs) << "," << target_string << "," << title_name << "," << to_string(dice_score[j]) << endl;
    }

}

void Read_Image_Size(){
    string buf,buf2, param;
	ifstream read_file(info_path);
	while(getline(read_file,buf)){

		if(buf.find("image size") != std::string::npos){
			istringstream iss(buf);
			iss >> param >> param >> param >> image_size[0] >> image_size[1] >> image_size[2];
            cout << "image size : " << image_size[0] << " " << image_size[1] << " " << image_size[2] << endl;
            logout << "image size : " << image_size[0] << " " << image_size[1] << " " << image_size[2] << endl;
		}
		if(buf.find("pixel dimension") != std::string::npos){
			istringstream iss(buf);
			iss >> param >> param >> param >>  pixel_spacing[0] >> pixel_spacing[1] >> pixel_spacing[2];
            cout << "pixel dimension : " << pixel_spacing[0] << " " << pixel_spacing[1] << " " << pixel_spacing[2] << endl;
            logout << "pixel dimension : " << pixel_spacing[0] << " " << pixel_spacing[1] << " " << pixel_spacing[2] << endl;

		}			
	}    
}


int main(int argc, char** argv){
    //cout << "version: " << CV_VERSION << std::endl;
    //getchar();
    stringstream log; //log用のss
    if(argc%2!= 1){
        cout << "invalid input" << endl;
        exit(1);
    }

    //引数から情報を受け取る
    for(int i = 1; i < argc; i+=2){
        string arg(argv[i]);
        string val(argv[i+1]);
        if(arg == "-outlog"){
            log << val;
            logout.open(log.str());
            cout << "Output log path : " << val << endl;
            logout << "Output log path : " << val << endl;   
        }
        else if(arg == "-iipath"){
            cout << "Import image path : " <<  val << endl;
            logout << "Import image path : " <<  val << endl;
            ii_path = val;
        }
        else if(arg == "-iiname"){
            cout << "Import image name : " << val << endl;
            logout << "Import image name : " << val << endl;
            ii_name = val;
        }
        else if(arg == "-ilpath"){
            cout << "Import label path : " <<  val << endl;
            logout << "Import label path : " <<  val << endl;
            il_path = val;
        }
        else if(arg == "-ilname"){
            cout << "Import label name : " << val << endl;
            logout << "Import label name : " << val << endl;
            il_name = val;
        }
        else if(arg == "-outimg"){
            cout << "Output tiff image path : " << val << endl;
            logout << "Output tiff image path : " << val << endl;
            out_img = val;
        }
        else if(arg == "-outscl"){
            cout << "Output scaling tiff image path : " << val << endl;
            logout << "Output scaling tiff image path : " << val << endl;
            out_scl = val;
        }
        else if(arg == "-outcsv"){
            cout << "Output csv path : " << val << endl;
            logout << "Output csv path : " << val << endl;
            out_csv = val;
        }
        else if(arg == "-trans"){
            trans_fp = val;
            cout << "Transform file path : " << val << endl;
            logout << "Transform file path : " << val << endl;
        }
        else if(arg == "-cs"){
            cs = stoi(val);
            cout << "Transform cross section : " << val << endl;
            logout << "Transform cross section : " << val << endl;
        }
        else if(arg == "-cs_x"){
            cs_x = stoi(val);
            cout << "Transform cross section x : " << val << endl;
            logout << "Transform cross section x : " << val << endl;
        }
        else if(arg == "-cs_y"){
            cs_y = stoi(val);
            cout << "Transform cross section y : " << val << endl;
            logout << "Transform cross section y : " << val << endl;
        }
        else if(arg == "-cs_z"){
            cs_z = stoi(val);
            cout << "Transform cross section z : " << val << endl;
            logout << "Transform cross section z : " << val << endl;
        }
        else if(arg == "-tar"){
            cout << "Target organ :";
            target_string = val;
            int tar_val = stoi(val);
            while(tar_val != 0){
                if(tar_val%10 == 1){
                    cout << " 1:Pancreas";
                    logout << " 1:Pancreas";
                    organ_flags[0] = true;
                }else if(tar_val %10 == 2){
                    cout << " 2:Stomach";
                    logout << " 2:Stomach";
                    organ_flags[1] = true;
                }else if(tar_val %10 == 3){
                    cout << " 3:Duodenum";
                    logout << " 3:Duodenum";
                    organ_flags[2] = true;
                }else if(tar_val %10 == 4){
                    cout << " 4:Left_Kidney";
                    logout << " 4:Left_Kidney";
                    organ_flags[3] = true;
                }

                tar_val /= 10;
            }
            cout << endl;
            logout << endl;
        }
        else if(arg == "-mode"){
            cout << "simulation mode :";
            logout << "simulation mode :";
            mode = stoi(val);
            if(mode == 1){
                cout << "Making data" << endl;
                logout << "Making data" << endl;
            }if(mode == 2){
                cout << "Estimation" << endl;
                logout << "Estimation" << endl;
            }
        }else if(arg == "-output"){
            if(stoi(val) == 1){
                cout << "in " << endl;
                cout << "image output mode" << endl;
                logout << "image output mode" << endl;
                OUTPUT_DATA_FLG = true;
            }
        }else if(arg == "-3Dtrans"){
            trans3D_fp = val;
            cout << "3D Transform file path : " << val << endl;
            logout << "3D Transform file path : " << val << endl;
        }else if(arg == "-out3Dcsv"){
            cout << "Output 3D csv path : " << val << endl;
            logout << "Output 3D csv path : " << val << endl;
            out3D_csv = val;
        }else if(arg == "-outElximg"){
            cout << "Output Elastix tiff image path : " << val << endl;
            logout << "Output Elastix tiff image path : " << val << endl;
            outElx_img = val;
        }
        else if(arg == "-aipath"){
            cout << "Answer image path : " <<  val << endl;
            logout << "Answer image path : " <<  val << endl;
            ai_path = val;
        }
        else if(arg == "-ainame"){
            cout << "Answer image name : " << val << endl;
            logout << "Answer image name : " << val << endl;
            ai_name = val;
        }
        else if(arg == "-dfile"){
            cout << "Output dice file : " << val << endl;
            logout << "Output dice file : " << val << endl;
            dice_file = val;
        }
        else if(arg == "-tarcsv"){
            cout << "Target csv path : " << val << endl;
            logout << "Target csv path : " << val << endl;
            target_csv = val;
        }else if(arg == "-info"){
            cout << "info path : " << val << endl;
            logout << "info path : " << val << endl;
            info_path = val;
        }else if(arg == "-unconverged"){
            cout << "unconverged list path : " << val << endl;
            logout << "unconverged list path : " << val << endl;
            unconverged_path = val;
        }else if(arg == "-th_val"){
            cout << "threshold val : " << val << endl;
            logout << "threshold val : " << val << endl;
            TH = stod(val);
        }else if(arg == "-divergence"){
            cout << "divergence list path : " << val << endl;
            logout << "divergence list path : " << val << endl;
            divergence_path = val;
        }else if(arg == "-initcsv"){
            cout << "initial 3D csv path : " << val << endl;
            logout << "initial 3D csv path : " << val << endl;
            init_csv = val;
        }else if(arg == "-df_csv"){
            cout << "deformation field csv: " << val << endl;
            logout << "deformation field csv: " << val << endl;
            df_csv = val;
        }else if(arg == "-cs_list"){
            cout << "cross section list: " << val << endl;
            logout << "cross section list: " << val << endl;
            cs_list = val;
        }else if(arg == "-df_ratio_list"){
            cout << "deformation ratio list: " << val << endl;
            logout << "deformation ratio list: " << val << endl;
            df_ratio_list = val;
        }else if(arg == "-VolumeAffine"){
            cout << "3D Affine file path : " << val << endl;
            logout << "3D Affine file path : " << val << endl;
            affine3D = val;
        }else if(arg == "-deformationField_path"){
            cout << "3D deformationField.raw path : " << val << endl;
            logout << "3D deformationField.raw path : " << val << endl;
            deformationField_3D = val;
        } else {
            cout << arg << " : Invalid input : " << val << endl;
            logout << arg << " : Invalid input : " << val << endl;
            exit(1);
        }
    }
    
    Read_Image_Size();
    InitGL();
    if(mode == 3) Read_Answer_Image();
    InitMPM();
    //object = new Particle_Cloud();//四角の箱が描かれるようになった（原さんはもともと書いていない）
    all_time_start = clock();

    //GLにより表示しない場合はここをコメントアウト
    
    while(1){
        Update();
    }
    //本当はここに　}
    //メインループ
    //描画(glutDisplayFunc ~ glutMainLoop) //InitGLも描画に必要 UpdateのglutPostRedisplayも
    /*
    glutDisplayFunc(Render);
    cout << "render :" << endl;
    glutIdleFunc(Update);
    glutMainLoop();
    */

    //delete fix;
    terminate_program();
    return 0;
}
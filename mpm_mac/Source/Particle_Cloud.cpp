#include "../Header/Particle_Cloud.h"
#include <iostream>
#include <omp.h>//</opt/homebrew/opt/libomp/include/omp.h>

Particle_Cloud::Particle_Cloud(){

}

Particle_Cloud::~Particle_Cloud(){

}

void Particle_Cloud::Draw(){
    glPointSize(4);
    glBegin(GL_POINTS);
    cout << "number of points: " << points.size() << endl;
    #pragma omp parallel for
    for(int i = 0; i < points.size(); i++){
        //膵臓：青
        //胃：赤
        //十二指腸：緑
        //腎臓：黄色
        if(points[i].dcmlabel == 1){
            glColor3f(0,0,1);
        }else if(points[i].dcmlabel == 2){
            glColor3f(1,0,0);
        }else if(points[i].dcmlabel == 3){
            glColor3f(0,1,0);
        }else if(points[i].dcmlabel == 4){
            glColor3f(1,1,0);
        }
        
        glVertex3d(points[i].position[0], points[i].position[1],points[i].position[2]);
        

        //駆動断面を黒で表す(目標位置)
        if(points[i].force_tag) {
            glColor3f(0,0,0);
            glVertex3d(points[i].target_position[0],points[i].target_position[1],points[i].target_position[2]);
        }


    }
    glEnd();

    GLenum err;
    while ((err = glGetError()) != GL_NO_ERROR) {
        std::cerr << "OpenGLエラー: " << err << std::endl;
    }

}

//応力を計算
void Particle_Cloud::Compute_Stress(){
    #pragma omp parallel for
    for(int p = 0; p < points.size(); p++){
        points[p].Compute_Stress();
        
    }
}

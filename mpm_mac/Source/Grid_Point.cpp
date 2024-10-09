#include "../Header/Grid_Point.h"

Grid_Point::Grid_Point(){
    mass = 0;
    force = velocity = momentum = Vector3d::Zero();
    //g = Vector3d(0.0,-9.8, 0.0);
    g = Vector3d::Zero();
    active = 0;
}

Grid_Point::~Grid_Point(){
    
}
void Grid_Point::Init(){

}
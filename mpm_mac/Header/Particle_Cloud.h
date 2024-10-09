#pragma once
#include "Particle.h"
#include <vector>
#include <GL/glut.h>

using namespace std;

class Particle_Cloud{
    public:
        Particle_Cloud();
        ~Particle_Cloud();

        void Compute_Stress();
        void Draw();
        void Init_Normal();

    public:
        vector<Particle> points;

};
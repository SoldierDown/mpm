#ifndef PARTICLE_H
#define PARTICLE_H
#include <eigen3/Eigen/Dense>

using namespace Eigen;

using Vec = Vector2d;
using Mat = Matrix2d;

class Particle
{
private:
    /* data */
public:
    Vec pos;
    Vec v;
    Vec p;
    double mass;
    double Jp;
    Mat F;
    Mat Fe;
    Mat C;
    Particle(Vec pos_, Vec v_ = Vec::Zero()): pos{pos_}, v{v_}, mass{1.0}, Jp{1.0},
                                                F{Mat::Identity()}, Fe{Mat::Identity()},
                                                C{Mat::Zero()}
    {
        p = mass * v;
    }
    ~Particle()
    {

    }
};



#endif
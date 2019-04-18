#ifndef PLANE_OBSTACLE_H
#define PLANE_OBSTACLE_H
#include <eigen3/Eigen/Dense>

using namespace Eigen;
using Vec = Vector2d;

class Plane_Obstacle
{
private:
    /* data */
public:
    Vec level;
    Vec normal;
    Vec tang;
    double mu;
    Plane_Obstacle(Vec level_, Vec normal_, Vec tang_): level{level_}, normal{normal_}, tang{tang_}, mu{0.5}{}
    ~Plane_Obstacle(){}
};



#endif
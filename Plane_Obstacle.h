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
    Plane_Obstacle(Vec level_, Vec normal_): level{level_}, normal{normal_}
    {
            
    }
    ~Plane_Obstacle(){}
};



#endif
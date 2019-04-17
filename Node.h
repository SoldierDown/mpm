#ifndef NODE_H
#define NODE_H
#include <eigen3/Eigen/Dense>

using namespace Eigen;
using Vec = Vector2d;

class Node
{
private:
    /* data */
public:
    double mass;
    Vec v;
    Vec p;
    Node(): mass{0.0}, v{Vec::Zero()}, p{Vec::Zero()}
    {

    }
    ~Node(){}
};


#endif
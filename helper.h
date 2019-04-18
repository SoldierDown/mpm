#ifndef HELPER_H
#define HELPER_H
#include <assert.h>
#include <eigen3/Eigen/Dense>
using namespace Eigen;
using Vec = Vector2d;
using Mat = Matrix3d;

double min(double a_, double b_)
{
    if(a_ < b_) return a_;
    else return b_;
}

double max(double a_, double b_)
{
    if(a_ > b_) return a_;
    else return b_;
}

double clamp(double x_, double a_, double b_)
{
    assert(a_ < b_);
    if(x_ < a_) return a_;
    else if(x_ > b_) return b_;
    else return x_;
}

Vec sqr(Vec vec_)
{
    return Vec(vec_(0) * vec_(0), vec_(1) * vec_(1));
}

Vec cub(Vec vec_)
{
    return Vec(pow(vec_(0), 3), pow(vec_(1), 3));
}

void polar_decomp(Mat F_, Mat& R_, Mat& S_)
{
    auto x = F_(0, 0) + F_(1, 1);
    auto y = F_(1, 0) - F_(0, 1);

    auto scale = 1.0 / sqrt(x * x + y * y);

    auto c = x * scale;        
    auto s = y * scale;

    R_(0, 0) = c;
    R_(0, 1) = -s;
    R_(1, 0) = s;
    R_(1, 1) = c;
    
    S_ = R_.transpose() * F_;    
}

#endif
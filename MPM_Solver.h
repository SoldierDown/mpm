#ifndef MPM_SOLVER_H
#define MPM_SOLVER_H
#include <vector>
#include <random>
#include "Node.h"
#include "Particle.h"
#include "Plane_Obstacle.h"
#include "helper.h"
#include "Constant.h"

using namespace std;

class MPM_Solver
{
private:

public:
    /* data */
    double E;
    double nu;
    double hardening;
    double dt;
    double frame_dt;

    double mu_0;
    double lambda_0;
    double dx;
    double inv_dx;
    double v_max;

    bool plastic;
    Vec gravity;

    vector<Node> nodes;
    vector<Particle> particles;
    vector<Plane_Obstacle> plane_obstacles;
    MPM_Solver(): E{1e4}, nu{0.2}, hardening{5}, dt{1e-4}, frame_dt{1e-3},
                            plastic{true}, gravity{Vec{0.0, -200.0}}                        
    {
        Init();
    }
    
    ~MPM_Solver(){}
    
    int Node_ID(int i_, int j_)
    {
        return i_ * (n + 1) + j_;
    }
    
    void Init()
    {
        Set_Parameters();
        Add_Objects(Vec(0.5, 0.5));
        Add_Obstacles();
        printf("Init()::SUCCESS!\n");
    }

    void Set_dt()
    {
        v_max = Get_Vmax();
        dt=max(1e-6, min(0.005, 0.1 * dx/max(v_max, 1e-2)));
    }

    double Get_Vmax()
    {
        double v_max = -1.0;
        for(auto &p: particles)
        {
            if(p.v.norm() > v_max) v_max = p.v.norm();
        }
        return v_max;
    }
    void Set_Parameters()
    {
        mu_0 = E / (2 * (1 + nu));
        lambda_0 = E * nu / ((1 + nu) * (1 - 2 * nu));
        dx = 1.0 / n;
        inv_dx = n / 1.0;
    }

    void Add_Obstacles()
    {
        Plane_Obstacle ground(Vec(0.0, 0.05), Vec(0.0, 1.0), Vec(1.0, 0.0));
        plane_obstacles.push_back(ground);
        Plane_Obstacle wall(Vec(0.05, 0.0), Vec(1.0, 0.0), Vec(0.0, 1.0));
        plane_obstacles.push_back(wall);
    }

    void Add_Objects(Vec center_)
    {
        double F;
        double r;
        double theta;
        for(int i = 0; i < NUMBER_OF_PARTICLE; i++)
        {
            F = rand()/(double)RAND_MAX;
            r = RADIUS * (double)sqrt(F);
            theta = TWO_PI * rand() / RAND_MAX;
            particles.push_back(Particle(center_ + Vec(r * cos(theta), r * sin(theta)), Vec(-20.0, 0.0)));
        }

        for(int i = 0; i <= n; i++)
        {
            for(int j = 0; j <= n; j++)
            {
                nodes.push_back(Node());
            }
        }
    }

    void Step()
    {
        //Set_dt();
        Reset_Node();
        Rasterize();
        Collision_Detection();
        Advection();
        printf("Step()::SUCCESS!\n");
    }

    void Reset_Node()
    {
        for(auto &node: nodes)
        {
            node.mass = 0.0;
            node.v.setZero();
            node.p.setZero();
        }
    }

    void Rasterize()
    {
        for(auto &p: particles)
        {
            Vector2i base_coord3 = (p.pos * inv_dx - Vec(1.0, 1.0)).cast<int>();
            Vec fx3 = p.pos * inv_dx - base_coord3.cast<double>();
            Vec w3[4] = {
                1.0 / 6.0 * cub(Vec(2.0, 2.0) - fx3),
                0.5 * cub(fx3 - Vec(1.0, 1.0)) - sqr(fx3 - Vec(1.0, 1.0)) + Vec(2.0/3.0, 2.0/3.0),
                0.5 * cub(Vec(2.0, 2.0) - fx3) - sqr(Vec(2.0, 2.0) - fx3) + Vec(2.0/3.0, 2.0/3.0),
                1.0 / 6.0 * cub(fx3 - Vec(1.0, 1.0))                
            };
            //printf("e: %f\n", exp(hardening * (1 - p.Jp)));
            auto e = max(HARDENING_MAX, exp(hardening * (1.0 - p.Jp)));
            auto mu = mu_0 * e;
            auto lambda = lambda_0 * e;

            auto J = p.F.determinant();

            Mat Re, Se;
            polar_decomp(p.Fe, Re, Se);

            auto Dinv = 3 * inv_dx * inv_dx;

            auto PF = 2 * mu * (p.Fe - Re) * p.Fe.transpose() + lambda * (J - 1) * J * Mat::Identity();

            auto vol = 1.0;
            auto stress = -(dt * vol) * (Dinv * PF);

            auto affine = stress + p.mass * p.C;

            int id;
            int ix, iy;
            for(int i = 0; i < 4; i++)
            {
                for(int j = 0; j < 4; j++)
                {
                    auto dpos = (Vec(i, j) - fx3) * dx;
                    ix = base_coord3.x() + i;
                    iy = base_coord3.y() + j;
                    id = Node_ID(ix, iy);

                    auto weight3 = w3[i].x() * w3[j].y();
                    nodes[id].mass += weight3 * p.mass;
                    //nodes[id].v += weight3 * (p.mass * p.v + affine * dpos);
                    nodes[id].p += weight3 * p.mass * p.v + weight3 * (affine * dpos);
                    
                    nodes[id].v = nodes[id].p / nodes[id].mass;
                }
            }

        }
    }

    void Collision_Detection()
    {
        for(int i = 0; i <= n; i++)
        {
            for(int j = 0; j <= n; j++)
            {
                auto &node = nodes[Node_ID(i, j)];
                if(node.mass > 0.0)
                {
                    node.p += dt * node.mass * gravity;
                    node.v = node.p / node.mass;
                    //node.v /= node.mass;
                    //node.mass = 1.0;
                    node.v += dt * gravity;
                    double x = (double) i / n;
                    double y = (double) j / n;

                    // double boundary = 0.05;
                    // if(x < boundary || x > 1 - boundary || y > 1 - boundary)
                    // {
                    //     node.mass = 0.0;
                    //     node.v.setZero();
                    //     node.p.setZero();
                    // }
                    // if(y < boundary)
                    // {
                    //     int sign = (node.v(0) > 0.0 ? 1: -1);
                    //     node.v(0) = sign * (max(0.0, abs(node.v(0)) - 1 * abs(node.v(1))));
                    //     node.v(1) = max(0.0, node.v(1));
                    //     node.p = node.mass * node.v;
                    // }
                    Vec pos(x, y);
                    for(auto &pl: plane_obstacles)
                    {
                        if(pl.normal.dot(pos - pl.level) < 0)
                        {
                            double vn = node.v.dot(pl.normal);
                            double vt = node.v.dot(pl.tang);
                            int sign = vt > 0.0 ? 1: -1;
                            double vnc = max(0.0, vn);
                            double vtc = sign * (max(0.0, abs(vt) - pl.mu * abs(vn)));
                            node.v = vnc * pl.normal + vtc * pl.tang;
                            if(vnc == 0.0 && vtc == 0.0) node.mass = 0.0;
                            node.p = node.mass * node.v;
                        }
                    }



                }
            }
        }
    }

    void Advection()
    {
        for(auto &p: particles)
        {
            Vector2i base_coord3 = (p.pos * inv_dx - Vec(1.0, 1.0)).cast<int>();
            Vec fx3 = p.pos * inv_dx - base_coord3.cast<double>();
            Vec w3[4] = {
                1.0 / 6.0 * cub(Vec(2.0, 2.0) - fx3),
                0.5 * cub(fx3 - Vec(1.0, 1.0)) - sqr(fx3 - Vec(1.0, 1.0)) + Vec(2.0/3.0, 2.0/3.0),
                0.5 * cub(Vec(2.0, 2.0) - fx3) - sqr(Vec(2.0, 2.0) - fx3) + Vec(2.0/3.0, 2.0/3.0),
                1.0 / 6.0 * cub(fx3 - Vec(1.0, 1.0))
            };    

            p.C.setZero();
            p.v.setZero();

            int ix, iy;
            int id;

            for(int i = 0; i < 4; i++)
            {
                for(int j = 0; j < 4; j++)
                {
                    auto dpos = (Vec(i, j) - fx3) * dx;
                    ix = base_coord3.x() + i;
                    iy = base_coord3.y() + j;
                    id = Node_ID(ix, iy);
                    auto node_v = nodes[id].v;
                    auto weight3 = w3[i].x() * w3[j].y();
                    p.v += weight3 * node_v;
                    p.C += 4 * inv_dx * ((weight3 * node_v) * dpos.transpose()); // HERE
                    
                }
            }


            p.pos += dt * p.v;
            Mat F = (Mat::Identity() + dt * p.C) * p.F;
            
            Mat svd_u, sig, svd_v;
            JacobiSVD<Mat> svd(F, ComputeFullU | ComputeFullV);
            svd_u = svd.matrixU();
            svd_v = svd.matrixV();
            sig =  svd_u.inverse() * F * svd_v.transpose().inverse();

            for(int i = 0; i < 2 * int(plastic); i++)
            {
                sig(i,i) = clamp(sig(i,i), 1.0 - 2.5e-2, 1.0 + 7.5e-3);
            }


            double oldJ = F.determinant();
            F = svd_u * sig * svd_v.transpose();
           
            //double Jp_new = clamp(p.Jp * oldJ / F.determinant(), 0.6, 20.0);
            double Jp_new = p.Jp * oldJ / F.determinant();
            p.Jp = Jp_new;
            p.F = F;       

        }
    }

};








#endif
#include "stdafx.h"

#include "point.hpp"

namespace geom {
point4d_t
  lerp( double lambda ,
        const point4d_t& p1,
        const point4d_t& p2)
  {
    point4d_t res;
    res.p[0] = (1-lambda) * p1.p[0] + lambda * p2.p[0];
    res.p[1] = (1-lambda) * p1.p[1] + lambda * p2.p[1];
    res.p[2] = (1-lambda) * p1.p[2] + lambda * p2.p[2];
    res.p[3] = (1-lambda) * p1.p[3] + lambda * p2.p[3];
    return res;
  }

  point3d_t
  lerp( double lambda ,
        const point3d_t& p1,
        const point3d_t& p2)
  {
    point3d_t res;
    res.p[0] = (1-lambda) * p1.p[0] + lambda * p2.p[0];
    res.p[1] = (1-lambda) * p1.p[1] + lambda * p2.p[1];
    res.p[2] = (1-lambda) * p1.p[2] + lambda * p2.p[2];
    return res;
  }

  point2d_t
  lerp( double lambda ,
        const point2d_t& p1,
        const point2d_t& p2)
  {
    point2d_t res;
    res.p[0] = (1-lambda) * p1.p[0] + lambda * p2.p[0];
    res.p[1] = (1-lambda) * p1.p[1] + lambda * p2.p[1];
    return res;
  }

  double
  lerp( double lambda ,
        const double& p1,//0
        const double& p2)//1
  {
    double res = (1-lambda) * p1 + lambda * p2;
    return res;
  }
  ////////////////////////////////
  point4d_t
  lerp( double lambda , double mu,
        const point4d_t& p1,//1,0
        const point4d_t& p2,//0,1
        const point4d_t& p3//0,0
        )
  {
    point4d_t res;
    res.p[0] = lambda * p1.p[0] + mu * p2.p[0] + (1 - lambda - mu) * p3.p[0];
    res.p[1] = lambda * p1.p[1] + mu * p2.p[1] + (1 - lambda - mu) * p3.p[1];
    res.p[2] = lambda * p1.p[2] + mu * p2.p[2] + (1 - lambda - mu) * p3.p[2];
    res.p[3] = lambda * p1.p[3] + mu * p2.p[3] + (1 - lambda - mu) * p3.p[3];
    return res;
  }

  point3d_t
  lerp( double lambda , double mu,
        const point3d_t& p1,
        const point3d_t& p2,
        const point3d_t& p3
        )
  {
    point3d_t res;
    res.p[0] = lambda * p1.p[0] + mu * p2.p[0] + (1 - lambda - mu) * p3.p[0];
    res.p[1] = lambda * p1.p[1] + mu * p2.p[1] + (1 - lambda - mu) * p3.p[1];
    res.p[2] = lambda * p1.p[2] + mu * p2.p[2] + (1 - lambda - mu) * p3.p[2];
    return res;
  }

  point2d_t
  lerp( double lambda , double mu,
        const point2d_t& p1,
        const point2d_t& p2,
        const point2d_t& p3
        )
  {
    point2d_t res;
    res.p[0] = lambda * p1.p[0] + mu * p2.p[0] + (1 - lambda - mu) * p3.p[0];
    res.p[1] = lambda * p1.p[1] + mu * p2.p[1] + (1 - lambda - mu) * p3.p[1];
    return res;
  }
  

}
//-*- mode:c++ -*-
#include "periodic_bspline.hpp"
#include "point.hpp"
template struct geom::periodic_bspline<geom::point2d_t>;
template struct geom::periodic_bspline<geom::point3d_t>;
template struct geom::periodic_bspline<geom::point4d_t>;
template struct geom::periodic_bspline<double>;

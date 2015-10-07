//-*- mode:c++ -*-
#include "periodic_bspline.hpp"
#include "point.hpp"
#include <numeric>
#include "tol.hpp"
#include "rmat.hpp"

template <class Point>
bool geom::periodic_bspline<Point>::check_invariants() const{

    auto pr = param_range();

    auto const & v1 = spl.eval_derivatives(degree()-1, pr.first );
    auto const & v2 = spl.eval_derivatives(degree()-1, pr.second );

    bool its_periodic = true;
    rmat<Point> rm(spl.control_points(), spl.knots(), degree() );
    double mult = rm.mult(knots()[degree()]);

    for(int i =0 ; i< degree() + 1 - mult; ++i){
        if(!tol::pt_eq(v1[i],v2[i] ))
        {
            its_periodic = false;
            break;
        }
    }

    return its_periodic;
}

template class geom::periodic_bspline<geom::point2d_t>;
template class geom::periodic_bspline<geom::point3d_t>;
template class geom::periodic_bspline<geom::point4d_t>;
template class geom::periodic_bspline<double>;

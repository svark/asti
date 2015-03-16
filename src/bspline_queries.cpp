#include "stdafx.h"
#include <numeric>
#include <vector>
#include "rmat.hpp"
#include "bspline_queries.hpp"

namespace geom {

//{{{ -- check is periodic
template <class SplineCurve>
bool
bspline_ops::is_periodic(const SplineCurve & crv)
{
    typedef typename SplineCurve::point_t point_t;
    enum {dim = point_dim<point_t>::dimension};

    int p = crv.degree();
    auto pb = crv.control_points().cbegin();
    auto pe = crv.control_points().cend();
    auto tb =  crv.knots().cbegin();
    auto te =  crv.knots().cend();

    size_t np = std::distance(pb, pe);
    size_t nt = std::distance(tb, te);
    if( np < 2)
        return false;

    if( nt != np + p + 1 )
        return false;

    if(nt < size_t(2*p + 2) )
        return false;

    auto toleq = [](const point_t &v,const point_t& w)  {
        return tol::eq(len(v-w),0);
    };

    if(!std::equal(pb, pb + p, pe - p, toleq))
        return false;

    std::vector<double> buf1(p + 1), buf2(p + 1);
    te -= (p + 1);

    auto deltas_at_start = buf1.begin();
    // fill up (1 + p) knot ranges) at the start
    std::adjacent_difference(tb,
                             tb + p + 1,
                             deltas_at_start);

    if( std::all_of( deltas_at_start + 1,
                     deltas_at_start + p + 1,
                     [](double v)->bool{
                         return tol::param_eq(v, 0);
                     }  ))
        return false;

    auto deltas_at_end = buf2.begin();
    // fill up 1 + p knot ranges at the end
    std::adjacent_difference(te,
                             te + p + 1,
                             deltas_at_end);

    if( std::all_of( deltas_at_end + 1,
                     deltas_at_end + p + 1,
                     [](double v)->bool{
                         return tol::param_eq(v, 0);
                     }  ))
        return false;

    return std::equal(deltas_at_start + 1,
                      deltas_at_start + p,
                      deltas_at_end + 1);
}
//}}}

//{{{ -- is bezier
template <class SplineType>
bool bspline_ops::is_bezier(const SplineType& spl)
{
    auto &t = spl.knots();
    int sz = spl.degree() + 1;
    if( t.size() != 2*sz || spl.control_points().size() != sz)
        return false;
    if(std::equal(t.begin(), t.begin() + sz,
                   t.rbegin() + t.size() - sz,  tol::param_eq ))
    {
        if(std::equal(t.rbegin(), t.rbegin() + sz,
                   t.begin() + t.size() - sz,  tol::param_eq ))
        {
            return true;
        }
    }
    return false;
}
//}}}

}

//{{{  instantiation scripts

/*
  Local Variables:
  eval:(load-file "./scripts/temp.el")
  eval:(setq methods (list "is_periodic" "is_bezier"
  ))
  eval:(setq spltypes (list "bspline<double>"
  "bspline<point2d_t>"
  "bspline<point3d_t>"
  "bspline<point4d_t>"
  ))
  eval:(instantiate-templates "bspline_queries" "bspline_ops" '()
  ((cons (car methods) spltypes ) ( cons (cadr methods) spltypes) ) )
  End:
// dump all explicitly instantiated templates below
*/
//}}}

//{{{  instantiation
#include "bspline.hpp"
#include "periodic_bspline.hpp"
#include "point.hpp"
namespace geom {
#include "bspline_queries_inst.cpp"
}
//}}}

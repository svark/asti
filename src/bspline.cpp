//--*-mode:c++-*-
//{{{ --includes
#include "stdafx.h"
#include <iostream>
#include <map>
#include <cmath>

#include "bspline.hpp"
#include "rmat.hpp"
#include "tol.hpp"
#include "util.hpp"
#include <utility>
#include <memory>
//}}}
template <class Point> struct bspline;

namespace geom{
//{{{ --(@* "evaluate curve at param @u")
template <class Point>
Point bspline<Point>::eval(double u) const
{
    if(tol::eq(u,t.back(),tol::param_tol))
        return cpts.back();

    rmat<Point> m(cpts,t, deg);
    return m.eval(u) + origin;
}
//}}}

//{{{ --(@* "evaluate derivative of the curve at the param @u")
// throws a spline_exception if u is not in parameter range of this curve
template <class Point>
typename bspline<Point>::vector_t
bspline<Point>::eval_derivative(int numDer, double u) const
{
    rmat<Point> m(cpts, t, deg);
    assert(numDer >= 0);
    return  make_vec( m.eval_derivative(numDer, u) ) ;
}

//}}}

//{{{ --(@* "parameter range")
template <class Point>
std::pair<double,double>
bspline<Point>::param_range() const
{
    return std::make_pair(t.front(), t.back() );
}
//}}}

//{{{ --(@* "constructors")

template <class Point>
bspline<Point>::bspline(cpts_t&& pts,
                        knots_t &&ks, int degree_):
    cpts(std::move(pts)),t(std::move(ks)),deg(degree_),
    origin(vector_t(0.0))
{
}

template <class Point>
bspline<Point>::bspline(std::tuple<cpts_t&&,knots_t&&,int&&> dat)
    :  cpts(std::forward<cpts_t>(std::get<0>(dat))),
       t(std::forward<knots_t>(std::get<1>(dat))),
       deg(std::get<2>(dat)) {}

template <class Point>
bspline<Point>::bspline(bspline<Point>&& other):
    cpts(std::forward<cpts_t>(other.cpts)),t(std::forward<knots_t>(other.t)),
    deg(other.deg),
    origin(std::forward<vector_t>(other.origin))
{
}

template <class Point>
void bspline<Point>::optimize()
{
    auto center = make_vec(centroid( cpts.begin(), cpts.end() ) );
    std::transform(cpts.begin(), cpts.end(), cpts.begin(),
                   [&center](Point& p){ return p -= center;} );
    origin += center;
}
//}}}

//{{{ --(@* "swap")
template <class Point>
void bspline < Point >::swap(bspline & other)
{
    t.swap(other.t);
    cpts.swap(other.cpts);
    std::swap(origin,other.origin);
    std::swap(deg, other.deg);
}
//}}}

//{{{ --(@*"evaluate bspline blossom f[0],f[1],...,f[deg-1]")
//see
// ./media/blossom1.png
// ./media/blossom2.png
template <class Point>
template <class KnotIter>
Point bspline<Point>::blossom_eval(KnotIter f)
{
    int p = degree();
    std::unique_ptr<point_t[]> cache(new point_t[p+1]);
    rmat<Point> m(cpts,t,deg);
    size_t nu = m.locate_nu(f[0]);
    for(int j = 0;j < p + 1; ++j) {
        cache[j] = control_points()[nu - p + j];
    }
    m.spline_compute(f, nu, cache.get());
    return cache[0] + origin;
}
//}}}
}
//{{{ (@* "template instantiations")

#include "point.hpp"


/*
  Local Variables:
  eval:(load-file "./temp.el")
  eval:(instantiate-templates "bspline" '("double"
  "point2d_t" "point3d_t"  "point4d_t")  '("blossom_eval") '("const double *"
  "std::vector<double>::const_iterator") )
  End:
*/
namespace geom {
#include "bspline_inst.cpp"
}
//}}}

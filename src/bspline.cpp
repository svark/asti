//--*-mode:c++-*-
//{{{ --includes
#include <iostream>
#include <map>
#include <cmath>

#include "bspline.hpp"
#include "rmat.hpp"
#include "tol.hpp"
#include "util.hpp"
#include <utility>
#include <memory>
#include "skip_idx_iter.hpp"
#include "transformed_back_inserter.hpp"
//}}}
namespace geom {
//{{{ --(@* "evaluate curve at param @u")
template <class Point>
Point bspline<Point>::eval(double u) const
{
    rmat<Point> m(cpts,t, deg);
    return m.eval(u) + origin;
}
//}}}

//{{{ --(@* "evaluate derivative of the curve at the param @u")

// throws a spline_exception if u is not in parameter range of this curve
template <class Point>
typename bspline<Point>::vector_t
bspline<Point>::eval_derivative(int derOrder, double u) const
{
    rmat<Point> m(cpts, t, deg);
    assert(derOrder >= 0);
    return make_vec( m.eval_derivative(derOrder, u) );
}

template <class Point>
typename bspline<Point>::vcpts_t
bspline<Point>::eval_derivatives(int derOrder, double u) const
{
    rmat<Point> m(cpts, t, deg);
    assert(derOrder >= 0);
    vcpts_t v;
    m.eval_derivatives(derOrder, u,
        util::transformed_back_inserter(v, &make_vec<point_t> ) );
    return v;
}
//}}}

//{{{ --(@* "parameter range")
template <class Point>
std::pair<double,double>
bspline<Point>::param_range() const
{
    size_t ncpts = control_points().size();
    return std::make_pair(t[deg], t[ncpts]);
}
//}}}

//{{{ --(@* "constructors")

template <class Point>
bspline<Point>::bspline(cpts_t  pts,
                        knots_t ks, int degree):
    t(std::move(ks)),cpts(std::move(pts)),deg(degree),
    origin(vector_t(0.0))
{
}

template <class Point>
bspline<Point>::bspline(std::tuple<cpts_t,knots_t,int>&& dat)
    :  t(std::forward<knots_t>(std::get<1>(dat))),
       cpts(std::forward<cpts_t>(std::get<0>(dat))),
       deg(std::get<2>(dat)),origin(vector_t(0.0))
{
}

template <class Point>
bspline<Point>::bspline(const bspline<Point>& other):
    t(other.t),
    cpts(other.cpts),
    deg(other.deg),
    origin(other.origin)
{
}

template <class Point>
bspline<Point>::bspline(bspline<Point>&& other):
    t(std::forward<knots_t>(other.t)),
    cpts(std::forward<cpts_t>(other.cpts)),
    deg(other.deg),
    origin(other.origin)
{
}
//}}}

//{{{  .(@* "optimize")
// optimize by storing all points relative to their centroid
template <class Point>
bspline<Point>& bspline<Point>::optimize()
{
    auto center = make_vec(centroid(cpts.cbegin(), cpts.cend()));
    std::transform(cpts.begin(), cpts.end(), cpts.begin(),
                   [&center](decltype(cpts[0]) p){ return p -= center;} );
    origin += center;
    return * this;
}

// deoptimize so that origin is at zero
template <class Point>
bspline<Point>& bspline<Point>::deoptimize()
{
    if(origin == vector_t(0.0))
        return *this;

    std::transform(cpts.begin(), cpts.end(), cpts.begin(),
                   [this](decltype(cpts[0]) p){ return p += this->origin;} );
    return * this;
}
//}}}

//{{{ --(@* "swap")
template <class Point>
void bspline < Point >::swap(bspline & other)
{
    if(& other == this)
        return;
    t.swap(other.t);
    cpts.swap(other.cpts);
    std::swap(origin,other.origin);
    std::swap(deg, other.deg);
}

//}}}

//{{{ --(@*"evaluate bspline blossom f[0],f[1],...,f[deg-1]")
//see
//(@file :file-name "./media/blossom1.png" :to "./media/blossom1.png" :display "blossom1")
//(@file :file-name "./media/blossom2.png" :to "./media/blossom2.png" :display "blossom2")
template <class Point>
template <class KnotIter>
Point bspline<Point>::blossom_eval(KnotIter f) const {
    int p = degree();
    std::unique_ptr<point_t[]> cache(new point_t[p+1]);
    rmat<Point> m(cpts,t,deg);
    size_t nu = m.locate_nu(*f);
    for(int j = 0;j < p + 1; ++j) {
        cache[j] = cpts[nu - p + j];
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
eval:(load-file "./scripts/temp.el")
eval:(instantiate-templates "bspline" "bspline" (list "double"
"point2d_t" "point3d_t"  "point4d_t")  (list (cons "blossom_eval" (list "const double *"
"std::vector<double>::const_iterator" "util::skip_ith_iter<std::vector<double>::const_iterator>" ) )))
End:
*/

namespace geom {
template <class Point> struct bspline;
#include "bspline_inst.inl"
}
//}}}

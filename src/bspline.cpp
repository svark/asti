//--*-mode:c++-*-
//{{{ --includes
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
    return rmat<Point>(cpts,t,deg).eval(u);
}
//}}}

//{{{ --(@* "evaluate derivative of the curve at the param @u")

template <class Point>
typename bspline<Point>::vector_t
bspline<Point>::eval_derivative(int derOrder, double u) const
{
    assert(derOrder >= 0);
    return make_vec(
        rmat<Point>(cpts, t, deg).eval_derivative(derOrder, u));
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
    t(std::move(ks)),cpts(std::move(pts)),deg(degree)
{
}

template <class Point>
bspline<Point>::bspline(std::tuple<cpts_t,knots_t,int>&& dat)
    :  t(std::forward<knots_t>(std::get<1>(dat))),
       cpts(std::forward<cpts_t>(std::get<0>(dat))),
       deg(std::get<2>(dat))
{
}

template <class Point>
bspline<Point>::bspline(const bspline<Point>& other):
    t(other.t),
    cpts(other.cpts),
    deg(other.deg)
{
}

template <class Point>
bspline<Point>::bspline(bspline<Point>&& other):
    t(std::forward<knots_t>(other.t)),
    cpts(std::forward<cpts_t>(other.cpts)),
    deg(other.deg)
{
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
    std::swap(deg, other.deg);
}

//}}}

//{{{ --(@*"evaluate bspline blossom at f[0],f[1],...,f[deg-1]")
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
    return cache[0];
}


template <class Point>
bool bspline<Point>::check_invariants() const
{
    if(deg <= 0)
        return false;

    if(t.size() != cpts.size()+deg+1) {
        return false;
    }

    if(cpts.size() < 2) {
        return false;
    }

    for(auto const &pt : cpts)
    {
        if(plen(pt) > 1 / tol::resabs)
            return false;
    }

    std::vector<double> uniqts(t);
    auto end = std::unique(uniqts.begin(),
                           uniqts.end(),
                           tol::param_eq);

    if( std::distance(uniqts.begin(),end)<=1 ) {
        return false;
    }

    return true;
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
//_

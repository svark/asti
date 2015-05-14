#include <numeric>
#include <vector>
#include "rmat.hpp"
#include "point_dim.hpp"
#include "constant_iterator.hpp"
#include "bspline_queries.hpp"
#include <limits>
#include "point.hpp"
namespace geom {
//{{{ -- check is periodic
template <class SplineCurve>
bool
ops::is_periodic(const SplineCurve & crv)
{
    typedef typename SplineCurve::point_t point_t;
    enum {dim = point_dim<point_t>::dimension};

    int p   = crv.degree();
    auto pb = crv.control_points().cbegin();
    auto pe = crv.control_points().cend();
    auto tb =  crv.knots().cbegin();
    auto te =  crv.knots().cend();

    size_t np = std::distance(pb, pe);
    size_t nt = std::distance(tb, te);

    if( np < 2 )
        return false;

    if( nt != np + p + 1 )
        return false;

    if(nt < size_t(2*p + 2) )
        return false;

    auto toleq = [](const point_t &v,const point_t& w)  {
        return tol::eq(len(v - w),0);
    };

    if(!std::equal(pb, pb + p, pe - p, toleq))
        return false;

    std::vector<double> buf1(p + 1), buf2(p + 1);
    te -= 2*(p + 1);
	tb += (p + 1);
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
bool ops::is_bezier(const SplineType& spl)
{
    auto &t = spl.knots();
    int sz = spl.degree() + 1;
    if( t.size() != 2*sz || spl.control_points().size() != sz)
        return false;
    if(std::equal(t.begin(), t.begin() + sz,
                  t.rbegin() + t.size() - sz, tol::param_eq ))
    {
        if(std::equal(t.rbegin(), t.rbegin() + sz,
                      t.begin() + t.size() - sz, tol::param_eq ))
        {
            return true;
        }
    }
    return false;
}
//}}}
//{{{ -- is regular

// determine if the first d + 1, and last d + 1 knots are equal
template <class SplineCurve>
bool
is_regular(const SplineCurve & c)
{
    auto & ts = c.knots();
    int p = c.degree();
    bool yes = std::equal(ts.cbegin(), ts.cbegin() + p + 1,
                           util::make_constant_iterator(ts[0]));
    yes &= std::equal(ts.crbegin(), ts.crbegin() + p + 1,
                      util::make_constant_iterator(ts[ts.size() - 1]));
    return yes;
}
//}}}
//{{{ -- curvature and torsion of a spline curve
template <class SplineType>
double ops::curvature(const SplineType & spl, double u)
{
    auto const & ds = spl.eval_derivatives(2, u);
    double t = len(ds[1]);
    if(tol::eq(t, 0))
		return std::numeric_limits<double>::infinity();
    return len(cross(ds[1], ds[2])) / (t * t * t);
}

template <class SplineType>
double ops::torsion(const SplineType & spl, double u)
{
    auto const &ds = spl.eval_derivatives(3, u);
    auto cp =  cross(ds[1], ds[2]);
    double w =  sqlen(cp);
	if(tol::eq(w, 0, tol::sqresabs))
        return std::numeric_limits<double>::infinity();
    double t = dot(decltype(cp)(ds[3]), cp);
    return t / w;
}

//}}}
}
//{{{  instantiation scripts

/*
  Local Variables:
  eval:(load-file "./scripts/temp.el")
  eval:(setq methods (list "is_periodic" "is_bezier" "curvature" "torsion"
  ))
  eval:(setq spltypes (list "bspline<double>"
  "bspline<point2d_t>"
  "bspline<point3d_t>"
  "bspline<point4d_t>"
  ))
  eval:(instantiate-templates "bspline_queries" "ops" (list )
          (product methods spltypes) )
  End:
// dump all explicitly instantiated templates below
*/
//}}}
//{{{  instantiation
#include "bspline.hpp"
#include "periodic_bspline.hpp"
#include "point.hpp"
namespace geom {
#include "bspline_queries_inst.inl"
}
//}}}

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
    auto pr = crv.param_range();
	
    auto const & v1 = crv.eval_derivatives(crv.degree()-1, pr.first );
    auto const & v2 = crv.eval_derivatives(crv.degree()-1, pr.second );
    bool its_periodic = true;
    for(int i =0 ; i< crv.degree(); ++i){
        if( !tol::eq(len(v1[i]-v2[i]), 0) )
        {
            its_periodic = false;
            break;
        }
    }
    return its_periodic;
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

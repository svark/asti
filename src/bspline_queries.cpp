#include <numeric>
#include <vector>
#include "rmat.hpp"
#include "point_dim.hpp"
#include "constant_iterator.hpp"
#include "bspline_queries.hpp"
#include <limits>
#include "periodic_bspline_cons.hpp"
#include "point.hpp"
namespace geom {
//{{{ -- check is periodic
template <class SplineCurve>
bool
qry::is_periodic(const SplineCurve & crv)
{
    auto const &c = get_spline(crv);
    return check_invariants(make_periodic_bspline(c.control_points(), c.knots(), c.degree()));
}

//}}}
//{{{ -- is bezier
template <class SplineType>
bool qry::is_bezier(const SplineType& spl)
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
//{{{ -- is clamped

// determine if the first d + 1, and last d + 1 knots are equal
template <class SplineCurve>
bool qry::is_clamped(const SplineCurve & c)
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
double qry::curvature(const SplineType & spl, double u)
{
    auto const & ds = spl.eval_derivatives(2, u);

    double t = len(ds[1]);
    if(tol::small(t))
        return std::numeric_limits<double>::quiet_NaN();

    return len(cross(ds[1], ds[2])) / (t * t * t);
}

template <class SplineType>
double qry::torsion(const SplineType & spl, double u)
{
    auto const & ds = spl.eval_derivatives(3, u);
    auto cp  =  cross(ds[1], ds[2]);
    double w =  sqlen(cp);
    if(tol::small(w, tol::sqresabs))
        return std::numeric_limits<double>::quiet_NaN();
    double t = dot(decltype(cp)(ds[3]), cp);
    return t / w;
}

//}}}
//{{{ auto lift dim
point3d_t qry::auto_lift_dim3(const point2d_t& p1, polynomial_tag, polynomial_tag)
{
    return point3d_t(p1);
}

point4d_t
qry::auto_lift_dim3(const point2d_t& p1, polynomial_tag, rational_tag)
{
    point4d_t p4d(p1);
    p4d[3] = 1.0;
    return p4d;
}

point4d_t
qry::auto_lift_dim3(const point3d_t& p1, polynomial_tag, rational_tag)
{
    return point4d_t(p1,1.0);
}


point4d_t
qry::auto_lift_dim3(const point3d_t& p1, rational_tag, rational_tag)
{
    point4d_t p4d(p1);
    std::swap(p4d[3],p4d[2]);
    return p4d;
}

point3d_t
qry::auto_lift_dim3(const point3d_t& p1, polynomial_tag, polynomial_tag)
{
    return p1;
}

point4d_t
qry::auto_lift_dim3(const point4d_t& p1, rational_tag, rational_tag)
{
    return p1;
}

point3d_t
qry::auto_lift_dim3(const point4d_t& p1, rational_tag, polynomial_tag)
{
    return scaled_copy(point3d_t(p1),1/p1[3]);
}

point3d_t
qry::auto_lift_dim3(const point3d_t& p1, rational_tag, polynomial_tag)
{
    return scaled_copy(point3d_t(point2d_t(p1)),1/p1[2]);
}


//}}}
}
//{{{  instantiation scripts

/*
  Local Variables:
  eval:(load-file "./scripts/temp.el")
  eval:(setq methods (list "is_periodic" "is_bezier" "is_clamped"
  "curvature" "torsion"
  ))
  eval:(setq spltypes (list "bspline<double>"
  "bspline<point2d_t>"
  "bspline<point3d_t>"
  "bspline<point4d_t>"
  ))
  eval:(instantiate-templates "bspline_queries" "qry" (list )
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

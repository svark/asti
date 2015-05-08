#include "smat.hpp"
#include "rmat.hpp"
#include "subdivide_curve.hpp"
#include "insert_knot.hpp"
#include "rotate_base_knot.hpp"
#include "trim_extend_join.hpp"
#include "rational_bspline_cons.hpp"
namespace geom {

template <class Point>
std::pair<bspline<Point>, bspline<Point>>
ops::split_open_curve(const bspline<Point>& spl, double u)
{
    return std::make_pair(clamp_at_right(u,spl),
                          clamp_at_left(u,spl));
}

template <class Point>
bspline<Point>
ops::split_periodic_curve(const periodic_bspline<Point>& pspl, double u)
{
    double s, e;
    int p = pspl.degree();
    std::tie(s, e) = pspl.param_range();
    if(tol::eq(pspl.periodic_param(u), s))
        return extract_regular_curve(pspl.spline());
    auto t = pspl.spline().knots();
    rmat_base_vd r(t, p);
    size_t nu = r.locate_nu(u);
    if(!tol::param_eq(t[nu] ,u))
        return split_periodic_curve(insert_knot(pspl, u),u);

    auto const & pc = rotate_base_knot(pspl, nu);
    return extract_regular_curve(pc.spline());
}



}


//{{{  instantiation scripts

/*
  Local Variables:
  eval:(load-file "./scripts/temp.el")
  eval:(setq methods (list "split_open_curve"
  "split_periodic_curve"
  ))
  eval:(setq spltypes (list "double"
  "point2d_t"
  "point3d_t"
  "point4d_t"
  ))
  eval:(instantiate-templates "subdivide_curve" "ops" (list )
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
#include "subdivide_curve_inst.inl"
}
//}}}

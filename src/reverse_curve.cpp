//-*- mode:c++ -*-
#include <vector>
#include "reverse_curve.hpp"
#include "bspline_x_cons.hpp"
#include "type_utils.hpp"
namespace geom{
//{{{ --(@* "reverse curve sense")

template <class SplineType>
SplineType ops::reverse_curve(SplineType spl)
{
    auto &cpts=spl.control_points();
    auto &t = spl.knots();
    size_t numPts = cpts.size();
    typedef typename SplineType::point_t point_t;
    typedef RAWTYPE(spl.control_points()) cpts_t;

    auto mod_fn = [] (typename cpts_t::iterator cb,
                      typename cpts_t::iterator ce,
                      typename std::vector<double>::iterator tb,
                      typename std::vector<double>::iterator te)
    {
        std::reverse(cb, ce);
        std::reverse(tb, te);
        std::transform(tb, te, tb,std::negate<double>());
    };

    typedef spline_traits<SplineType> str;
    return SplineType(std::move(spl), mod_fn);
}

//}}}

}

//{{{  instantiation scripts

/*
  Local Variables:
  eval:(load-file "./scripts/temp.el")
  eval:(setq methods (list "reverse_curve"
  ))
  eval:(setq spltypes (list "bspline<double>"
  "bspline<point2d_t>"
  "bspline<point3d_t>"
  "bspline<point4d_t>"
  "periodic_bspline<double>"
  "periodic_bspline<point2d_t>"
  "periodic_bspline<point3d_t>"
  "periodic_bspline<point4d_t>"
  "rational_bspline < point2d_t,regular_tag>"
  "rational_bspline < point3d_t,regular_tag>"
  "rational_bspline < double, regular_tag>"
  "rational_bspline < point2d_t,periodic_tag>"
  "rational_bspline < point3d_t,periodic_tag>"
  "rational_bspline < double,periodic_tag>"
  ))
  eval:(instantiate-templates "reverse_curve" "ops" (list )
  (product methods spltypes) )
  End:
  // dump all explicitly instantiated templates below
  */
//}}}

//{{{  instantiation
#include "bspline.hpp"
#include "periodic_bspline.hpp"
#include "rational_bspline.hpp"
#include "point.hpp"
namespace geom {
#include "reverse_curve_inst.inl"
}
//}}}

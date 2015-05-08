#include <vector>
#include "reverse_curve.hpp"
#include "bspline_x_cons.hpp"
#include "type_utils.hpp"
namespace geom{
//{{{ --(@* "reverse curve sense")

template <class SplineType>
SplineType ops::reverse_curve(const SplineType& spl)
{
    auto &cpts=spl.control_points();
    auto &t = spl.knots();
    size_t numPts = cpts.size();
    typedef typename SplineType::point_t point_t;
    typedef RAWTYPE(spl.control_points()) cpts_t;
    cpts_t new_cpts(numPts);

    std::vector<double> new_knots(spl.knots().size());

    std::reverse_copy(cpts.cbegin(), cpts.cend(), new_cpts.begin());
    std::reverse_copy(t.cbegin(), t.cend(), new_knots.begin());

    std::transform(new_knots.begin(), new_knots.end(),
                   new_knots.begin(),
                   std::negate<double>() );

    return make_bsplinex<SplineType>(std::move(new_cpts),
                                     std::move(new_knots),
                                     spl.degree()
        ).translate(spl.base_point());
}

template <class SplineType>
SplineType& ops::inplace_reverse_curve(SplineType& spl)
{
    std::reverse(spl.control_points().begin(), spl.control_points().end());
    std::reverse(spl.knots().begin(), spl.knots().end());

    std::transform(spl.knots().begin(), spl.knots()().end(),
                   spl.knots().begin(),
                   std::negate());
    return spl;
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
  "rational_bspline<bspline<point2d_t>>"
  "rational_bspline<bspline<point3d_t>>"
  "rational_bspline<bspline<point4d_t>>"
  "rational_bspline<periodic_bspline<point2d_t>>"
  "rational_bspline<periodic_bspline<point3d_t>>"
  "rational_bspline<periodic_bspline<point4d_t>>"
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
#include "point.hpp"
namespace geom {
#include "reverse_curve_inst.inl"
}
//}}}

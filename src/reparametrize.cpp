//-*- mode:c++ -*-
#include "point.hpp"
#include "reparametrize.hpp"
#include "tol.hpp"
#include "geom_exception.hpp"
#include "bspline_x_cons.hpp"
#include "spline_traits.hpp"
namespace geom {

//{{{ --(@* "re-parametrise curve with a new range")
template <class SplineType>
SplineType ops::reparametrize(const SplineType& spl,
                              double t1, double t2)
{
    typedef typename SplineType::knots_t knots_t;
    typedef typename SplineType::point_t point_t;
    knots_t new_knots( spl.knots().size() );
    auto & t = spl.knots();

    double first_t,last_t ;
    std::tie(first_t,last_t)= spl.param_range();

    if( fabs(last_t - first_t) < tol::param_tol)
        throw geom_exception(bad_knot_spacing_t);

    const double scale = 1.0 / ( last_t - first_t);
    size_t i = 0;
    for(double u : t )
    {
        const double par = ( u - first_t) * scale;
        new_knots[i++] = t1 + (t2 - t1) * par;
    }
    typedef spline_traits<SplineType> str;

    typedef RAWTYPE(spl.control_points()) cpts_t;
    return make_bspline(std::move(cpts_t(spl.control_points())),
                        std::move(new_knots),
                        spl.degree(),
                        typename str::ptag(),
                        typename str::rtag());
}

template <class SplineType>
SplineType ops::reparametrize_start(const SplineType& spl,
                                    double t1)
{

    double first_t,last_t ;
    std::tie(first_t,last_t)= spl.param_range();
    double t2 = t1 + (last_t - first_t);
    return reparametrize(spl, t1, t2);
}
//}}}

}

/*
  Local Variables:
  eval:(load-file "./scripts/temp.el")
  eval:(setq methods (list "reparametrize" "reparametrize_start"
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
  eval:(instantiate-templates "reparametrize" "ops" (list )
  (product methods spltypes ))
  End:
  // dump all explicitly instantiated templates below
  */
//}}}

//{{{  instantiation
#include "bspline.hpp"
#include "periodic_bspline.hpp"
#include "point.hpp"
namespace geom {
#include "reparametrize_inst.inl"
}
//}}}

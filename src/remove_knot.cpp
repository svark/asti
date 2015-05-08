#include "bspline_x_cons.hpp"
#include "remove_knot.hpp"
#include "raise_degree.hpp"
#include "rmat.hpp"
#include "type_utils.hpp"
namespace geom {

template <typename SplineCurve>
SplineCurve
ops::fair_by_knot_removal(const SplineCurve & crv_, double tol)
{
    SplineCurve crv(crv_);
    while(crv.degree() < 3)
        crv.swap(raise_degree(crv));
    auto const &  ts =  crv.knots();
    std::vector<double> uts;
    uts.reserve(ts.size());
    std::unique_copy(ts.cbegin(),
                     ts.cend(),
                     uts.begin());
    int j = 0;
    std::vector<double> kdash_variation(uts.size());
    for(auto t : uts)
    {
        auto kdash_minus = crv.eval_derivative(3, t - tol::param_tol/ 2);
        auto kdash_plus  = crv.eval_derivative(3, t + tol::param_tol/ 2);
        kdash_variation[j] = sqlen(kdash_plus - kdash_minus);
        ++j;
    }

    auto it = std::max_element(kdash_variation.cbegin(),
                               kdash_variation.cend());
    size_t k = std::distance(kdash_variation.cbegin(), it);
    double u = uts[k];
    rmat_base_vd rm(crv.knots(), crv.degree());
    size_t nu = rm.locate_nu(u);
    auto const & cpts = crv.control_points();
    auto const & t = crv.knots();
    SplineCurve::cpts_t newcpts(cpts);
    typedef RAWTYPE(cpts[0]) point_t;
    auto l =  (t[nu + 1] - t[nu - 3]) * make_vec(cpts[nu - 1])  -
        (t[nu + 1] - t[nu]) * make_vec( cpts[nu - 2])  ;

    l *= 1.0 / (t[nu] - t[nu - 3]);

    auto r =  (t[nu + 3] - t[nu - 1]) * make_vec( cpts[nu + 1])   -
        (t[nu] - t[nu - 1]) * make_vec ( cpts[nu + 2] ) ;

    r *= 1.0 / (t[nu] - t[nu - 3]);

    newcpts[nu] = make_pt((t[nu + 2] - t[nu]) * l + (t[nu] - t[nu - 2]) * r);
    newcpts[nu] *= 1 / (t[nu + 2] - t[nu - 2]);

    double den = sqlen(newcpts[nu] - cpts[nu]);
    if(den > tol * tol) {
        newcpts[nu] = cpts[nu] + tol *
            (newcpts[nu] - cpts[nu]) / sqrt(den);
    }else {
        newcpts[nu] = cpts[nu];
    }

    return make_bsplinex<SplineCurve>(std::move(newcpts),
                         std::move(t), crv.degree() );
}

template <typename SplineCurve>
rational_bspline < SplineCurve >
ops::fair_by_knot_removal(const rational_bspline<SplineCurve> & crv,
                     double tol)
{
    return make_rbspline( fair_by_knot_removal(crv.spline(), tol));
}


}

/*
  Local Variables:
  eval:(load-file "./scripts/temp.el")
  eval:(setq methods (list "fair_by_knot_removal"
  ))
  eval:(setq spltypes (list "bspline<double>"
  "bspline<point2d_t>"
  "bspline<point3d_t>"
  "bspline<point4d_t>"
  "periodic_bspline<double>"
  "periodic_bspline<point2d_t>"
  "periodic_bspline<point3d_t>"
  "periodic_bspline<point4d_t>"
  "rational_bspline < bspline<point2d_t>>"
  "rational_bspline < bspline<point3d_t>>"
  "rational_bspline < bspline<point4d_t>>"
  "rational_bspline < periodic_bspline<point2d_t>>"
  "rational_bspline < periodic_bspline<point3d_t>>"
  "rational_bspline < periodic_bspline<point4d_t>>"
  ))
  eval:(instantiate-templates "remove_knot" "ops" (list ) (product
  methods spltypes ))
  End:
// dump all explicitly instantiated templates below
*/
#include "bspline.hpp"
#include "periodic_bspline.hpp"
#include "rational_bspline.hpp"
#include "point.hpp"
namespace geom {
#include "remove_knot_inst.inl"
}

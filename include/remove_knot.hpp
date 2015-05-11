#ifndef ASTI_REMOVE_KNOT
#define ASTI_REMOVE_KNOT
#include "tol.hpp"
#include "bspline_fwd.hpp"
namespace geom { namespace ops {
// farin cagd page 422
template <typename SplineCurve>
extern SplineCurve fair_by_knot_removal(const SplineCurve & crv,
                                        double tol);

template <typename Point,class PTag>
rational_bspline<Point, PTag >
fair_by_knot_removal(const rational_bspline < Point, PTag > & crv,
                     double tol)
{
  return make_rbspline( fair_by_knot_removal(crv.spline(), tol));
}
}
}
#endif // ASTI_REMOVE_KNOT

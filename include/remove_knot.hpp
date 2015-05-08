#ifndef ASTI_REMOVE_KNOT
#define ASTI_REMOVE_KNOT
#include "tol.hpp"
#include "bspline_fwd.hpp"
namespace geom { namespace ops {
// farin cagd page 422
template <typename SplineCurve>
extern SplineCurve fair_by_knot_removal(const SplineCurve & crv,
                                        double tol);

template <typename SplineCurve>
extern rational_bspline < SplineCurve >
fair_by_knot_removal(const rational_bspline < SplineCurve > & crv,
                     double tol);
}
}
#endif // ASTI_REMOVE_KNOT

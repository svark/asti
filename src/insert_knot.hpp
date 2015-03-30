#ifndef INSERT_KNOT
#define INSERT_KNOT
#include <vector>
#include "rational_bspline_cons.hpp"
namespace geom{
namespace ops {

template <class SplineType>
extern  SplineType
insert_knots(const SplineType& crv,
             const std::vector<double>& refined_knots);

template <class SplineType>
extern  SplineType
insert_knot(const SplineType& crv, double u);

template <class SplineCurve>
rational_bspline <SplineCurve>
insert_knot(const rational_bspline<SplineCurve>
            & crv, double u)
{
    return make_rbspline(insert_knot(crv.spline(), u));
}
}
}
#endif // INSERT_KNOT

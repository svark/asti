#ifndef RAISE_DEGREE
#define RAISE_DEGREE
#include "bspline_fwd.hpp"

namespace geom {
namespace ops
{

template <class SplineType>
extern SplineType
raise_degree(const SplineType&crv);

template <class SplineCurve>
extern rational_bspline <SplineCurve>
raise_degree(const rational_bspline<SplineCurve>
             & crv, double u);

}
}
#endif // RAISE_DEGREE

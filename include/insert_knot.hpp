#ifndef ASTI_INSERT_KNOT
#define ASTI_INSERT_KNOT
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
bool
match_knots(SplineCurve & spl1, SplineCurve & spl2);

}
}
#endif // ASTI_INSERT_KNOT

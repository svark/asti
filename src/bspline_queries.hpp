#include "bspline_fwd.hpp"
namespace geom {
namespace ops {

template <class SplineType>
extern bool is_bezier(const SplineType& spl);

template <class SplineType>
bool is_bezier(const rational_bspline < SplineType > & spl)
{
    return is_bezier(spl.spline());
}

template <class SplineType>
extern bool is_periodic(const SplineType & spl);

template <class SplineType>
bool is_periodic(const rational_bspline < SplineType >& spl)
{
    return is_periodic(spl.spline());
}
}}

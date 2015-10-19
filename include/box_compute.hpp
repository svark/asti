#ifndef ASTI_BOX_COMPUTE
#define ASTI_BOX_COMPUTE
#include "box.hpp"
#include "bspline_fwd.hpp"
namespace geom {
template <class P> class circle;
template <class P> class conic_arc;

namespace ops
{
template <class SplineType>
extern box<typename SplineType::point_t>
compute_box(const SplineType &spl);

template <class SplineType>
extern box<typename SplineType::point_t>
compute_box_tight(const SplineType &spl);

template <class Point>
extern box<Point>
compute_box(const circle<Point> &spl);

template <class Point>
extern box<Point>
compute_box(const conic_arc<Point> &c);

}
}
#endif // ASTI_BOX_COMPUTE

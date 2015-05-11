#ifndef ASTI_BOX_COMPUTE
#define ASTI_BOX_COMPUTE
#include "box.hpp"
#include "circle.hpp"
namespace geom {
    namespace ops
    {
        template <class SplineType>
        extern box<typename SplineType::point_t>
        compute_box(const SplineType &spl);

	template <class Point>
        extern box<Point>
        compute_box(const circle<Point> &spl);
    }
}
#endif // ASTI_BOX_COMPUTE

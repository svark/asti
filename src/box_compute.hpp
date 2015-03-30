#pragma once
#include "box.hpp"
namespace geom {
    namespace ops
    {
        template <class SplineType>
        extern box<typename SplineType::point_t>
        compute_box(const SplineType &spl);
    }
}

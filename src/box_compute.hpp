#pragma once
#include "box.hpp"
#include "spline_traits.hpp"
namespace geom { 
    namespace bspline_ops
    {
        template <class SplineType>
        extern box<spline_traits<SplineType>::dim>
        compute_box(const SplineType &spl);
    }
}
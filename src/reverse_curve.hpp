#pragma once
namespace geom {
namespace bspline_ops {
template <class SplineType>
    extern SplineType reverse_curve(SplineType& spl);

    template <class SplineType>
    extern SplineType& inplace_reverse_curve(SplineType& spl);
}
}
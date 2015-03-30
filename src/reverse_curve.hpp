#pragma once
namespace geom {
namespace ops {
template <class SplineType>
    extern SplineType reverse_curve(const SplineType& spl);

    template <class SplineType>
    extern SplineType& inplace_reverse_curve(SplineType& spl);
}
}

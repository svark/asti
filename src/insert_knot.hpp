#pragma once
#include <vector>
namespace geom{
    namespace bspline_ops {
        template <class SplineType>
        extern  SplineType
            insert_knots(const SplineType& crv,
            const std::vector<double>& refined_knots);

        template <class SplineType>
        extern  SplineType
            insert_knot(const SplineType& crv,double u);
        
    }
}
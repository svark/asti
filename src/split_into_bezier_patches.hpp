#pragma once
#include <list>
#include "Eigen/Core"
namespace geom {
namespace bspline_ops{
template <class SplineType>
    extern std::list<SplineType,Eigen::aligned_allocator<SplineType> >
    split_into_bezier_patches(const SplineType &spl);
   
}
}
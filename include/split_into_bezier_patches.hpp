//-*- mode: c++ -*-
#ifndef ASTI_SPLIT_INTO_BEZIER_PATCHES
#define ASTI_SPLIT_INTO_BEZIER_PATCHES

#include <list>
#include "Eigen/Core"
#include "bspline_fwd.hpp"
namespace geom {
namespace ops{
template <class SplineType>
extern std::list<SplineType,Eigen::aligned_allocator<SplineType> >
split_into_bezier_patches(const SplineType &spl);

template <class SplineType>
extern SplineType
first_bezier_patch(const SplineType &spl);

template <class SplineType>
extern SplineType
last_bezier_patch(const SplineType &spl);

}
}
#endif // ASTI_SPLIT_INTO_BEZIER_PATCHES

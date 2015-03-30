#ifndef SPLIT_INTO_BEZIER_PATCHES
#define SPLIT_INTO_BEZIER_PATCHES

#include <list>
#include "Eigen/Core"
#include "bspline_fwd.hpp"
namespace geom {
namespace ops{
template <class SplineType>
    extern std::list<SplineType,Eigen::aligned_allocator<SplineType> >
    split_into_bezier_patches(const SplineType &spl);

template <class SplineType>
extern std::list< rational_bspline < SplineType>,
           Eigen::aligned_allocator<rational_bspline < SplineType >>>
split_into_bezier_patches(const rational_bspline < SplineType > &spl);

}
}
#endif // SPLIT_INTO_BEZIER_PATCHES

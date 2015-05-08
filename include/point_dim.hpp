#ifndef ASTI_POINT_DIM
#define ASTI_POINT_DIM
#include <type_traits>
#include "Eigen/Core"
namespace geom{

template <class Point>
struct point_dim
{
    enum {dimension = Point::dimension};
    typedef Eigen::aligned_allocator<Point> alloc_t;
};

template <>
struct point_dim<double>
{
    enum {dimension = 1};
    typedef std::allocator<double> alloc_t;
};

template<class PointVec>
struct inc_dimension
{
    typedef decltype(higher_dim(PointVec())) type;
};

template<class PointVec>
struct dec_dimension
{
    typedef decltype(lower_dim(PointVec())) type;
};

}
#endif // ASTI_POINT_DIM

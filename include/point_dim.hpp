#ifndef ASTI_POINT_DIM
#define ASTI_POINT_DIM
#include <type_traits>
#include "point_fwd.hpp"
namespace geom{

template <class Point>
struct point_dim
{
};

template <int dim>
struct point_dim<pt_t < dim >>
{
    enum {dimension = dim};
    //    typedef Eigen::aligned_allocator<pt_t < dim > > alloc_t;
};

template <>
struct point_dim<double>
{
    enum {dimension = 1};
    //    typedef std::allocator<double> alloc_t;
};

template <int dim>
struct point_dim<vec_t < dim >>
{
    enum {dimension = dim};
    //typedef Eigen::aligned_allocator<vec_t < dim > > alloc_t;
};


template<class PointVec>
struct inc_dimension
{
    typedef decltype(higher_dim(std::declval<PointVec>())) type;
};

template<class PointVec>
struct dec_dimension
{
    typedef decltype(lower_dim(std::declval<PointVec>())) type;
};

}
#endif // ASTI_POINT_DIM

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
};

template <>
struct point_dim<double>
{
    enum {dimension = 1};
};

template <int dim>
struct point_dim<vec_t < dim >>
{
    enum {dimension = dim};
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

#define ENABLE_IF_DIM_IS_2(PointU) typename std::enable_if<point_dim < PointU >::dimension == 2, int>::type = 0
}
#endif // ASTI_POINT_DIM

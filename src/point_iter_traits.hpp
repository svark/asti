#ifndef POINT_ITER_TRAITS
#define POINT_ITER_TRAITS
//#include <type_traits>
#include "Eigen/Dense"
//#include "vtens.hpp"
#include <vector>
#include "point_dim.hpp"
namespace geom {
template <class PointIter>
struct point_iter_traits
{
    typedef typename std::iterator_traits<PointIter>::value_type point_t;
    enum {dim = point_dim<point_t>::dimension};

     typedef std::vector<point_t,
                         typename point_dim<point_t>::alloc_t> PointContT;
     typedef std::vector<decltype(make_vec(point_t())),
                         typename point_dim<point_t>::alloc_t>
                         VectorContT;
    //typedef typename std::conditional<dim==1,std::vector<double>,vector_of_type< point_t > >::type PointContT;
    //typedef typename std::conditional<dim==1,std::vector<double>,vector_of_type< decltype(make_vec(point_t())) > >::type VectorContT;

};



}
#endif // POINT_ITER_TRAITS

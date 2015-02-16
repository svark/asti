#ifndef SPLINE_TRAITS_HPP
#define SPLINE_TRAITS_HPP
#include "point.hpp"
//#include <boost/mpl/integral_c.hpp>
namespace geom {
typedef std::false_type regular_tag;
typedef std::true_type periodic_tag;
template <class Point> struct periodic_bspline;

template <class SplineType>
struct spline_traits_base
{
    typedef typename SplineType::point_t point_t;
    enum{dim= point_traits < point_t *>::dim};
};

template <class SplineType>
struct spline_traits: spline_traits_base < SplineType >
{
    typedef regular_tag ptag;
};
template <class Point>
struct spline_traits<periodic_bspline < Point >> : spline_traits_base < periodic_bspline < Point > >
{
    typedef periodic_bspline ptag;
};
}
#endif

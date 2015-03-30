#ifndef SPLINE_TRAITS_HPP
#define SPLINE_TRAITS_HPP
#include "point_iter_traits.hpp"
#include "bspline_fwd.hpp"

namespace geom {
typedef std::false_type regular_tag;
typedef std::true_type periodic_tag;

typedef std::false_type nonrational_tag;
typedef std::true_type rational_tag;

template <class SplineType>
struct spline_traits_base
{
    typedef typename SplineType::point_t point_t;
    enum{dim= point_iter_traits < point_t *>::dim};
    typedef SplineType type;
};

template <class SplineType>
struct spline_traits: spline_traits_base < SplineType >
{
    typedef regular_tag ptag;
    typedef nonrational_tag rtag;
};
template <class Point>
struct spline_traits<periodic_bspline < Point >>
    : spline_traits_base < periodic_bspline < Point > >
{
    typedef periodic_tag ptag;
    typedef nonrational_tag rtag;
};

template <class Point>
struct spline_traits<rational_bspline<periodic_bspline < Point > > >
    : spline_traits_base < periodic_bspline < Point > >
{
    typedef periodic_tag ptag;
    typedef rational_tag rtag;
};
template <class Point>
struct spline_traits<rational_bspline<bspline < Point > > >
    : spline_traits_base <rational_bspline< bspline < Point > > >
{
    typedef regular_tag ptag;
    typedef rational_tag rtag;
};


}
#endif

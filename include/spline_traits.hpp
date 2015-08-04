#ifndef ASTI_SPLINE_TRAITS_HPP
#define ASTI_SPLINE_TRAITS_HPP
#include "point_dim.hpp"
#include "bspline_fwd.hpp"

namespace geom {

struct regular_tag : std::false_type {};
struct periodic_tag : std::true_type {};

struct  polynomial_tag : std::false_type{};
struct  rational_tag: std::true_type{};

template <class SplineType>
struct spline_traits;

template<class PTag, class RTag, class Point>
struct get_traits_type_from_tags
{
};

template<class Point>
struct get_traits_type_from_tags < periodic_tag, polynomial_tag, Point >
{
    typedef spline_traits <  periodic_bspline < Point>> type;
};

template<class Point>
struct get_traits_type_from_tags < regular_tag, polynomial_tag,Point >
{
    typedef spline_traits <  bspline < Point > > type;
};

template<class Point>
struct get_traits_type_from_tags < periodic_tag, rational_tag, Point >
{
    typedef spline_traits <
        rational_bspline < Point, periodic_tag > > type;
};

template<class Point>
struct get_traits_type_from_tags<regular_tag, rational_tag, Point>
{
    typedef spline_traits < rational_bspline <  Point, regular_tag> > type;
};


template <class SplineType>
struct spline_traits_base
{
    typedef typename SplineType::point_t point_t;
    enum{dimension = point_dim < point_t >::dimension};
    typedef SplineType spline_type;
};

template <class SplineType>
struct spline_traits: spline_traits_base < SplineType >
{
    typedef regular_tag ptag;
    typedef polynomial_tag rtag;
};
template <class Point>
struct spline_traits<periodic_bspline < Point >>
    : spline_traits_base < periodic_bspline < Point > >
{
    typedef periodic_tag ptag;
    typedef polynomial_tag rtag;
};

template <class Point>
struct spline_traits<rational_bspline< Point,periodic_tag > >
    : spline_traits_base < periodic_bspline < Point > >
{
    typedef periodic_tag ptag;
    typedef rational_tag rtag;
};
template <class Point>
struct spline_traits<rational_bspline< Point, regular_tag > >
    : spline_traits_base <bspline<  Point > >
{
    typedef regular_tag ptag;
    typedef rational_tag rtag;
};

}
#endif

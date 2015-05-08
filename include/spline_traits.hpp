#ifndef ASTI_SPLINE_TRAITS_HPP
#define ASTI_SPLINE_TRAITS_HPP
#include "point_dim.hpp"
#include "bspline_fwd.hpp"

namespace geom {
typedef std::false_type regular_tag;
typedef std::true_type periodic_tag;

typedef std::false_type polynomial_tag;
typedef std::true_type rational_tag;

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
        rational_bspline < periodic_bspline < Point > >> type;
};

template<regular_tag, rational_tag, class Point>
struct get_traits_type_from_tags
{
    typedef spline_traits < rational_bspline < bspline < Point >>> type;
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

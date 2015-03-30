#ifndef SUBDIVIDE_CURVE_HPP
#define SUBDIVIDE_CURVE_HPP
#pragma once
#include "bspline.hpp"
#include "periodic_bspline.hpp"
#include "rational_bspline_cons.hpp"
#include <utility>//for pair
namespace geom {
namespace ops {

template <class Point>
extern std::pair<bspline<Point>,bspline<Point> >
split_open_curve(const bspline<Point>& spl, double u);

template <class Point>
extern bspline<Point>
split_periodic_curve(const periodic_bspline<Point>& pspl, double u);

template <class Point>
rational_bspline<bspline<Point> >
split_periodic_curve(const rational_bspline<periodic_bspline<Point> > & pspl, double u)
{
   return make_rbspline( split_periodic_curve( pspl.spline(), u ) );
}

template <class Point>
std::pair<rational_bspline<bspline<Point> >, rational_bspline<bspline<Point> > >
split_open_curve(const rational_bspline<bspline<Point> > & pspl, double u)
{
   return make_rbspline( split_open_curve( pspl.spline(), u ) );
}

}
}
#endif // SUBDIVIDE_CURVE_HPP

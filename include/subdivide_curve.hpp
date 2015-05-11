#ifndef ASTI_SUBDIVIDE_CURVE_HPP
#define ASTI_SUBDIVIDE_CURVE_HPP
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
rational_bspline<Point, regular_tag >
split_periodic_curve(const rational_bspline<Point, periodic_tag > & pspl, double u)
{
   return make_rbspline( split_periodic_curve( pspl.spline(), u ) );
}

template <class Point>
std::pair<rational_bspline<Point, regular_tag >, rational_bspline<Point, regular_tag> >
split_open_curve(const rational_bspline<Point, regular_tag > & pspl, double u)
{
  auto& pr = split_open_curve( pspl.spline(), u );
  return std::make_pair(make_rbspline(pr.first),make_rbspline(pr.second));
}

}
}
#endif // ASTI_SUBDIVIDE_CURVE_HPP

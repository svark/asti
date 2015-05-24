#pragma once
#include "bspline.hpp"
namespace geom {
namespace ops {

template <class Fn>
extern bspline<double>
cubic_approx1d(Fn f, std::vector<double> t);

template <class FnType>
extern bspline<double>
quad_approx1d(FnType f, std::vector<double> t);

template <class SplineType>
extern double
foot_param(const SplineType &spl,
           const typename SplineType::point_t& p);

}
}

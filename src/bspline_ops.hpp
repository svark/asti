//-*- mode: c++ -*-
#ifndef _BSPLINE_OPS_H
#define _BSPLINE_OPS_H
#include "box.hpp"
#include "spline_traits.hpp"
namespace geom {
struct bspline_ops {
    template <class SplineType>
    static SplineType
    insert_knots(const SplineType& crv,
                 const std::vector<double>& refined_knots);

    template <class SplineType>
    static SplineType
    insert_knot(const SplineType& crv,
                double u);

    template <class SplineType>
    static SplineType
    raise_degree(const SplineType&crv);

    template <class SplineType>
    static std::list<SplineType>
    split_into_bezier_patches(const SplineType &spl);

    template <class SplineType>
    static std::list<SplineType>
    split_into_bezier_patches_hard(const SplineType &spl);

    template <class SplineType>
    static SplineType
    extract_curve(const SplineType &spl, double a, double b);

    template <class SplineType>
    static SplineType reverse_curve(SplineType& spl);

    template <class SplineType>
    static SplineType& inplace_reverse_curve(SplineType& spl);

    template <class SplineType>
    static SplineType reparametrize(const SplineType& spl,
                                    double t1, double t2);

    template <class SplineType>
    static bool is_bezier(const SplineType& spl);

    template <class PointIter, class KnotIter>
    static bool
    is_periodic(PointIter, PointIter,
                KnotIter, KnotIter, int degree);

    template <class Fn>
    static bspline<double>
    cubic_approx1d(Fn f, std::vector<double>& t);

    template <class FnType>
    static bspline<double>
    quad_approx1d(FnType f, std::vector<double>& t);

    template <class SplineType>
    static double
    foot_param(const SplineType &spl,
               const typename SplineType::point_t& p);

    template <class SplineType>
    static box<spline_traits<SplineType>::dim>
    compute_box(const SplineType &spl);
};
}
#endif //_BSPLINE_OPS_H

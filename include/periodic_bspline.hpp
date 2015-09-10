#ifndef ASTI_PERIODIC_BSPLINE_HPP
#define ASTI_PERIODIC_BSPLINE_HPP

#include "bspline.hpp"
#include <utility>
#include <vector>
#include "util.hpp"

namespace geom
{

static
double periodic_param(const std::pair<double, double>& range,
                      double u)
{
    if (u <_in_> range)
        return u;

    double s = range.first, e = range.second;
    double r = fmod(u - s, (e - s));
    if( r <  0 ) r += (e - s);
    return s + r;
}


template <class Point>
struct periodic_bspline
{
    typedef Point   point_t;
    typedef decltype(Point() - Point())  vector_t;
    enum {dimension = point_dim<point_t>::dimension};
    typedef typename bspline<Point>::cpts_t cpts_t;
    typedef typename bspline<Point>::knots_t knots_t;
    typedef typename bspline<Point>::vcpts_t vcpts_t;
    typedef std::tuple<cpts_t,knots_t, int> tuple_t;

    periodic_bspline(cpts_t pts, knots_t ks, int degree_)
        :spl(std::move(pts),std::move(ks),degree_)
    {
    }

    periodic_bspline(bspline<point_t>&& other)
        :spl(std::forward<bspline<point_t>>(other))
    {
    }

    periodic_bspline(const bspline<point_t>& other)
        :spl(other)
    {
    }

    periodic_bspline(periodic_bspline&& other)
        :spl(std::forward<bspline<point_t> >(other.spl))
    {
    }

    template <class ModifierFn>
    periodic_bspline(periodic_bspline&& other, ModifierFn mod_fn)
        :spl(std::forward<bspline<point_t> >(other.spl),mod_fn)
    {
    }

    periodic_bspline(const periodic_bspline& other)
        :spl(other.spl)
    {
    }

    double periodic_param(double u) const
    {
        return geom::periodic_param(param_range(),u);
    }

    point_t eval(double u) const {
        return spl.eval(periodic_param(u));
    }

    const bspline<Point>& spline() const { return spl; }

    template <class KnotIter>
    point_t blossom_eval(KnotIter f) const {
        std::vector<double> modfs(degree());
        double d = periodic_param(*f) - *f;

        std::transform(f,
                       f + degree(),
                       modfs.begin(),
                       [&](double v) { return v+d; });
        return spl.blossom_eval(modfs.cbegin());
    }

    auto eval_derivative(int derOrder,double u) const -> RAWTYPE(std::declval<bspline<point_t>>().eval_derivative(derOrder,u))
    {
        return spl.eval_derivative(derOrder,periodic_param(u));
    }
    auto eval_derivatives(int derOrder,double u) const -> RAWTYPE(std::declval<bspline<point_t>>().eval_derivatives(derOrder,u))
    {
        return spl.eval_derivatives(derOrder,periodic_param(u));
    }

    std::pair<double,double> param_range() const {
        return spl.param_range();
    }

    void swap( periodic_bspline & other ) {
        spl.swap(other.spl);
    }

    // accessors
    const knots_t & knots() const { return spl.knots(); }
    const cpts_t &  control_points() const { return spl.control_points(); }
    int degree() const { return spl.degree(); };
private:
    bspline<point_t> spl;
};

template <class Point>
bool check_invariants(const periodic_bspline<Point> &spl) {
    return check_invariants(spl.spline());
}

}
#endif

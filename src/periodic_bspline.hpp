#ifndef PERIODIC_BSPLINE_HPP
#define PERIODIC_BSPLINE_HPP
#include "bspline.hpp"
#include <utility>
#include <vector>
#include "util.hpp"
namespace geom
{
struct spline_wrap_t{};

template <class SplineType>
double periodic_param(const SplineType & spl, double u)
{
    auto p  = spl.param_range();
    if (u <_in_> p)
        return u;

    double s = p.first, e = p.second;
    double r = fmod(u - s, (e - s));
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

    template <class PointU>
    static periodic_bspline < PointU >    rebound_type();

    /*static std::tuple<cpts_t,knots_t, int>
    wrap(const cpts_t& pts, const knots_t &ks, int degree);*/


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

    periodic_bspline(const periodic_bspline& other)
        :spl(other.spl)
    {
    }

    double periodic_param(double u) const
    {
        return geom::periodic_param(*this,u);
    }

    point_t eval(double u) const {
        return spl.eval(periodic_param(u));
    }

    const bspline<Point>& spline() const { return spl; }

    template <class knot_iter>
    point_t blossom_eval(knot_iter f) const {
        std::vector<double> modfs(degree());
        double d = periodic_param(*f) - *f;

        std::transform(f,
                       f + degree(),
                       modfs.begin(),
                       [&](double v) { return v+d; });
        return spl.blossom_eval(modfs.cbegin());
    }

    vector_t eval_derivative(int numDer,double u) const
    {
        return spl.eval_derivative(numDer,periodic_param(u));
    }
    vcpts_t eval_derivatives(int numDer,double u) const
    {
        return spl.eval_derivatives(numDer,periodic_param(u));
    }

    std::pair<double,double> param_range() const {
        return spl.param_range();
    }

    periodic_bspline&
    translate(const vector_t& t) {
        spl.translate(t); return *this;
    }

    /*constexpr*/ bool is_periodic() const { return true; }

    // store cpts relative to cg
    void optimize() { spl.optimize(); }

    void swap( periodic_bspline & other ) {
        spl.swap(other.spl);
    }

    // accessors
    const knots_t & knots() const { return spl.knots(); }
    const cpts_t &  control_points() const { return spl.control_points(); }
    int degree() const { return spl.degree(); };
    const vector_t& base_point() const { return spl.base_point(); }
private:
    bspline<point_t> spl;
};




}
#endif

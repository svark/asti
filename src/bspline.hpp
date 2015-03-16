#ifndef BSPLINE_H
#define BSPLINE_H
#include <algorithm>
#include <iterator>
#include <vector>
#include <utility>
#include <vector>
#include <list>
#include <Eigen/Core>
#include "point_iter_traits.hpp"
namespace geom {

template <class Point>
struct bspline {

    typedef Point   point_t;
    typedef decltype(Point() - Point())  vector_t;
    enum {dimension = point_dim<point_t>::dimension};

    template <class PointU>
    static bspline < PointU >    rebound_type();

    typedef typename point_iter_traits<point_t*>::PointContT  cpts_t;
    typedef typename point_iter_traits<point_t*>::VectorContT vcpts_t;
    typedef std::vector<double> knots_t;

    bspline(std::tuple<cpts_t,knots_t,int>&& dat);
    bspline(cpts_t pts, knots_t ks, int degree_);

    bspline(const bspline& other);
    bspline(bspline&& other);

    point_t eval(double u) const;

    template <class knot_iter>
    point_t blossom_eval(knot_iter f) const;

    vector_t eval_derivative(int numDer,double u) const;

    /// derivatives from 0->numDer at the param @u
    vcpts_t eval_derivatives(int numDer, double u) const;

    std::pair<double,double> param_range() const;

    bspline<Point>& translate(const vector_t& t) {
        origin += t; return *this;
    }

    // store cpts relative to cg
    void optimize();

    void swap( bspline & other ) ;

    // accessors
    const knots_t & knots() const { return t;}
    const cpts_t &  control_points() const { return cpts;}
    int degree() const { return deg; };
    const vector_t& base_point() const { return origin; }
protected:
    knots_t t;
    cpts_t cpts;
    int deg;
    vector_t origin;
};

template <class Point>
bool operator==(const bspline<Point>& bs1,const bspline<Point>&bs2) {
    return std::make_tuple(bs1.control_points(),bs1.knots(),bs1.base_point(),bs1.degree())
        ==std::make_tuple(bs2.control_points(),bs2.knots(),bs2.base_point(),bs2.degree());
}

}
#endif //BSPLINE_H

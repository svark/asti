#ifndef ASTI_BSPLINE_H
#define ASTI_BSPLINE_H
#include <algorithm>
#include <iterator>
#include <vector>
#include <utility>
#include <vector>
#include <memory>
//#include <Eigen/Core>
#include "point_dim.hpp"
#include "type_utils.hpp"
namespace geom {

template <class Point>
struct bspline {

    typedef Point   point_t;
    typedef decltype(make_vec(point_t()))  vector_t;
    enum {dimension = point_dim<point_t>::dimension};

    typedef decltype(mk_stdvec(point_t()))  cpts_t;
    typedef decltype(mk_stdvec(vector_t())) vcpts_t;

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

    void swap( bspline & other ) ;

    // accessors
    const knots_t&  knots()              const { return t;}
    const cpts_t &  control_points()     const {
        return cpts;
    }
    int             degree()             const { return deg; };

protected:

    knots_t t;
    cpts_t cpts;
    int deg;
};

template <class Point>
bool operator==(const bspline<Point>& bs1, const bspline<Point>& bs2) {

    bs1.swap(
        reparametrize(extract_regular_curve(bs1), 0, 1).deoptimize());
    bs2.swap(
        reparametrize(extract_regular_curve(bs2), 0, 1).deoptimize());

    if(!match_degrees(bs1, bs2))
        return false;
    if(!match_knots(bs1, bs2))
        return false;
    return bs1.control_points() == bs2.control_points();
}

}
#endif //BSPLINE_H

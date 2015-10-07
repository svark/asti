#ifndef ASTI_BSPLINE_H
#define ASTI_BSPLINE_H
#include <algorithm>
#include <iterator>
#include <vector>
#include <utility>
#include <vector>
#include <memory>
#include "point_dim.hpp"
#include "type_utils.hpp"
#include "tol.hpp"

namespace geom {

template <class Point>
class bspline {
public:
    typedef Point   point_t;
    typedef VECTOR_TYPE(Point)  vector_t;
    enum {dimension = point_dim<point_t>::dimension};

    typedef ARRAY_TYPE(point_t)  cpts_t;
    typedef ARRAY_TYPE(vector_t) vcpts_t;

    typedef std::vector<double> knots_t;

    bspline(std::tuple<cpts_t,knots_t,int>&& dat);
    bspline(cpts_t pts, knots_t ks, int degree_);

    bspline(const bspline& other);
    bspline(bspline&& other);

    template <class ModifierFn>
    bspline(bspline&& other, ModifierFn mod_fn)
    {
        bspline(std::forward<bspline>(other)).swap(*this);
        mod_fn(cpts.begin(),cpts.end(),t.begin(),t.end());
        check_invariants();
    }

    point_t eval(double u) const;

    template <class KnotIter>
    point_t blossom_eval(KnotIter f) const;

    vector_t eval_derivative(int numDer,double u) const;

    /// derivatives from 0->numDer at the param @u
    vcpts_t eval_derivatives(int numDer, double u) const;

    std::pair<double,double> param_range() const;

    void swap( bspline & other ) ;

    bool check_invariants() const;

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

    if(bs1.degree() != bs2.degree())
        return false;

    if(!std::equal(bs1.knots().cbegin(),bs1.knots().cend(),
                   bs2.knots().cbegin() , tol::param_eq) )
        return false;

    if(!std::equal(bs1.control_points().cbegin(),bs1.control_points().cend(),
                   bs2.control_points().cbegin() , tol::pt_eq<Point> ))
        return false;
    return true;
}

template <class Point>
bool check_invariants(const bspline<Point>& spl) {return spl.check_invariants(); }

}
#endif //BSPLINE_H

#ifndef _BSPLINE_H
#define _BSPLINE_H
#include <algorithm>
#include <iterator>
#include <vector>
#include <utility>
#include <vector>
#include <list>
#include <Eigen/Core>

namespace geom {

    template <class Point>
    struct bspline {

        typedef Point   point_t;
        typedef decltype(Point() - Point())  vector_t;

        typedef std::vector<point_t> cpts_t;
        typedef std::vector<double> knots_t;

        bspline(const cpts_t& pts, const knots_t &ks, int degree_);

        bspline(std::tuple<cpts_t&&,knots_t&&,int&&> dat)
            :  cpts(std::forward<cpts_t>(std::get<0>(dat))),
               t(std::forward<knots_t>(std::get<1>(dat))),
               deg(std::get<2>(dat)) {}

        bspline(cpts_t&& pts, knots_t &&ks, int degree_);

        bspline(const bspline& other);

        bspline(bspline&& other);

        point_t eval(double u) const;

        template <class knot_iter>
        point_t blossom_eval(knot_iter f);

        vector_t eval_derivative(int numDer,double u) const;

        std::pair<double,double> param_range() const;

        bspline<Point>& translate(const vector_t& t) {
            origin += t; return *this;
        }

        //constexpr
        static int dimension() {
            return point_traits<point_t*>::dim;
        }
        // store cpts relative to cg
        void optimize();

        void swap( bspline & other ) {
            t.swap(other.t);
            cpts.swap(other.cpts);
            std::swap(origin,other.origin);
        }

        // accessors
        const knots_t & knots() const { return t;}
        const cpts_t &  control_points() const { return cpts;}
        int degree() const { return deg; };
        const vector_t& base_point() const { return origin; }
    protected:
        knots_t t;
        cpts_t cpts;
        vector_t origin;
        int deg;
    };
	
	
}

#endif //_BSPLINE_H

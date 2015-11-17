#ifndef ASTI_CIRCLE_HPP
#define ASTI_CIRCLE_HPP
#include "point_dim.hpp"
#include <math.h> //sin and cos
#include "bspline_fwd.hpp"
#include "spline_traits.hpp"
#include "tol.hpp"
#include "type_utils.hpp"


namespace geom{
template <class Point>
class circle
{
public:
    enum {dim = point_dim<Point>::dimension};
    typedef Point point_t;
    typedef VECTOR_TYPE(point_t) vector_t;

    template <class PointU>
    circle(const PointU& center_,
           const PointU& point_,
           VECTOR_TYPE(PointU) const & ydir_)
        :origin(center_), px(point_), y(ydir_)
    {
        auto x = px - origin;
        assert(tol::not_small(len(x)));
        y -= x*dot(x,y)/len(x);
        y = normalize(y) * len(x);
        assert(check_invariants());
    }


    template <class PointU>
    circle(const PointU& center_, const PointU & point_,
           ENABLE_IF_DIM_IS_2(PointU))
        :origin(center_), px(point_)
    {
        auto x = px - origin;
        assert(tol::not_small(len(x)));
        y   = normalize( make_vec( - x[1], x[0]) ) * len(x);
    }


    point_t eval(double u) const
    {
        auto x = (px - origin);
        return center() + x * cos(u) + y * sin(u);
    }

    vector_t tangent(double u) const
    {
        auto x = (px - origin);
        return  - x * sin(u) + y * cos(u);
    }

    vector_t normal(double u) const
    {
        return eval(u) - center();
    }
    // accessors

    decltype(cross(vector_t(),vector_t()))
    plane_normal() const
    {
        auto const &x = (px - origin);
        return  normalize( cross(x,y) );
    }

    vector_t ydir() const
    {
        return y;
    }

    double radius() const
    {
        return len(px - center());
    }

    point_t start() const { return px; }

    point_t center() const { return origin; }

    bool check_invariants() {
        if(tol::pt_eq(px,origin)) return false;
        return tol::small(dot(px - origin , y)) ;
    }
private:
    point_t  origin;
    point_t  px;
    vector_t y;//cached y direction
};


template <class Point>
extern double
foot_param(const circle<Point> &c,
           const Point& p);

template <class Point>
extern rational_bspline< Point, regular_tag >
make_rbspline_from_circle(const circle<Point>& circ);


template <class Point>
extern circle<Point>
make_circle(const Point& p1,
            const Point& p2,
            const Point& p3);
}




#endif

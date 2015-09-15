#ifndef ASTI_CIRCLE_HPP
#define ASTI_CIRCLE_HPP
#include "point_dim.hpp"
#include <math.h>
#include "bspline_fwd.hpp"
#include "spline_traits.hpp"
#include "tol.hpp"
#include "type_utils.hpp"


namespace geom{
template <class Point>
struct circle
{
    enum {dim = point_dim<Point>::dimension};
    typedef Point point_t;
    typedef decltype(make_vec(point_t())) vector_t;

    template <class PointU>
    circle(const PointU& center_,
           const PointU& point_,
           VECTOR_TYPE(PointU) const & ydir_)
        :center(center_), start_pt(point_), ydir(ydir_)
    {
        auto x = start_pt - center;
        assert(tol::not_small(len(x)));
        ydir -= x*dot((start_pt-center),ydir)/len(x);
        ydir = normalize(ydir) * len(x);
    }


    template <class PointU>
    circle(const PointU& center_, const PointU & point_,
           ENABLE_IF_DIM_IS_2(PointU))
        :center(center_), start_pt(point_)
    {
        auto x = start_pt - center;
        ydir   = normalize( make_vec( - x[1], x[0]) ) * len(x);
    }


    //  (@file :file-name "media/circle2.png" :to "./media/circle2.png" :display "eval at param")
    point_t eval(double u) const
    {
        auto x = (start_pt - center);
        auto y = ydir;
        return center + x * cos(u) + y * sin(u);
    }

    vector_t tangent(double u) const
    {
        auto x = (start_pt - center);
        auto y = ydir;
        return  - x * sin(u) + y * cos(u);
    }

    vector_t normal(double u) const
    {
        return eval(u) - center;
    }
    // accessors
    point_t  getCenter() const
    {
        return center;
    }

    point_t getStart() const
    {
        return start_pt;
    }

    decltype(cross(vector_t(),vector_t()))
    getPlaneNormal() const
    {
        auto const &x = (start_pt - center);
        auto const &y = ydir;
        return  normalize( cross(x,y) );
    }

    vector_t getYDir() const
    {
        return ydir;
    }

    double getRadius() const
    {
        return len(start_pt - center);
    }
private:
    point_t  center;
    point_t  start_pt;
    vector_t ydir;
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

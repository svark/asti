#ifndef ASTI_LINE_HPP
#define ASTI_LINE_HPP

#include <algorithm>
#include "point_fwd.hpp"
#include "geom_exception.hpp"
#include "tol.hpp"
namespace geom{

template <class Point>
struct line
{
    typedef Point point_t;
    typedef decltype(make_vec(Point())) vector_t;

    line(const point_t& p1,
         const vector_t& v)
        :start(p1), dir(normalize(v))
    {
        assert(tol::eq(len(dir) , 1.0));
    }

    point_t eval(double u) const
    {
        return start + dir * u;
    }
    point_t start_pt() const { return start; }
    vector_t direction() const { return dir ;};
private:
    point_t  start;
    vector_t dir;
};

template <class Point>
static line<Point>
make_line(const Point& p1,
          const Point& p2
    )
{
    auto dir = normalize(p2-p1);
    return line<Point>(p1, dir);
}

template <class Point>
struct line_seg
{
    typedef Point point_t;
    typedef decltype(Point() - Point()) vector_t;

    line_seg(const point_t& p1,
             const point_t& p2)
        :l(p1,normalize(p2-p1)),a(0),b(len(p2-p1))
    {
    }

    line_seg(line<point_t> l_,
             double a_,
             double b_)
        :l(l_),a(a_),b(b_)
    {
    }

    // evaluate line at the given parameter
    point_t eval(double u) const { return l.eval(u); }

    // return the parameter range in which this segment lies
    std::pair<double,double> param_range() const {
        return std::make_pair(a,b);
    }

    line<point_t> getLine() const { return l;}

private:
    line<point_t> l;
    double a, b;
};

template <class Point>
static line_seg<Point>
make_line_seg(const Point& p1,
              const Point& p2
    )
{
    return line_seg<Point>(p1, p2);
}

//.  ./media/intersect_lines1.png

template <class Point>    
std::pair<Point, Point >
closest_points(const line < Point >& l1,
               const line < Point >& l2)
{
    auto d1 = l1.direction();
    auto d2 = l2.direction();
    auto r =  l1.start_pt() - l2.start_pt();
    double b = dot(d1, d2);
    double c = dot(d1, r);
    double f = dot(d2, r);

    double d = 1.0 - b * b;

    if(tol::eq(d, 0))
        return std::make_pair(l1.start_pt(), l2.start_pt()+f*l2.direction());
    
    auto p1 = l1.start_pt() +  d1 * (dot(d1, (d2 * f - r )) / d);
    auto p2 = l2.start_pt() +  d2 * (dot(d2, (r - d1 * c)) / d);
    return std::make_pair(p1, p2);
}

extern point2d_t
intersect_lines(const line < point2d_t >& l1,
                const line < point2d_t >& l2);

}

#endif

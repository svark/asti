//-*- mode:c++ -*-
#include "circle.hpp"
#include "geom_exception.hpp"
#include "rational_bspline_cons.hpp"
#include "smat.hpp"
// using geom::method convention to enable the script below to
// instantiate

using geom::circle;
using geom::regular_tag;
using geom::rational_bspline;

template <class Point>
double
geom::foot_param(const circle<Point> &c,
           const Point& p)
{
    auto x = normalize(c.getStart() - c.getCenter());
    auto y = normalize(c.getYDir());

    auto v = (p - c.getCenter());

    if(tol::small(dot(v,y)) && tol::small(dot(v,x)))
        throw geom_exception(point_at_axis_error);

    return atan2(dot(v,y) , dot(v,x) );
}

template <class Point>
circle<Point>
geom::make_circle(const Point& p1,
                  const Point& p2,
                  const Point& p3)
{
    typedef Point point_t;
    auto v1 = p1 - p2;
    const double eps = tol::resabs;

    double lv1 = sqlen(v1);
    if(lv1 < eps*eps)
        throw geom_exception(circle_too_small);

    auto v2 = p2 - p3;
    double lv2 = sqlen(v2);
    if(lv2 < eps*eps)
        throw geom_exception(circle_too_small);

    auto v3 = p3 - p1;
    double lv3 = sqlen(v3);

    double d = sqlen( cross( v1, v2) );
    auto center = dot(v1, -v3) / d*d;
    double lens[] = {lv1, lv2, lv3};
    decltype(v1) vs[] = {v1,v2,v3};

    size_t idx = std::max_element(lens,lens+3) - lens;

    v1 = vs[(idx + 1)%3];   lv1 = lens[(idx + 1)%3];
    v2 = vs[(idx + 2)%3];   lv2 = lens[(idx + 2)%3];
    v3 = vs[idx];           lv3 = lens[idx];

    point_t pts[] = {p1,p2,p3};
    std::rotate(pts, pts + (idx + 1)%3, pts + 3);
    auto normal  = cross(v1,v2);
    double denom = sqlen(normal);
    if(denom < eps*eps )
        throw geom_exception(circle_too_small);

    double alpha = lv2 * dot(v1,-v3) *0.5 / ( denom );
    double beta  = lv3 * dot(-v1,v2) *0.5 / ( denom );

    auto cp = lerp(alpha,beta,pts[0],pts[1],pts[2]);
    decltype(normal) xdir(pts[0] - cp);
    auto nnormal (normalize(normal));
    assert(!isnan(nnormal));
    return circle<point_t>(cp, pts[0], decltype(v1)(cross(nnormal, xdir)));
}

template <class Point>
rational_bspline< Point, regular_tag >
geom::make_rbspline_from_circle(const circle<Point>& circ)
{
    auto start_pt =  circ.getStart();
    auto center =  circ.getCenter();
    auto const & x = (start_pt - center);
    decltype(x)  y = cross(circ.getPlaneNormal(), decltype(circ.getPlaneNormal())(x) );

    double radius =  circ.getRadius();

    auto a = start_pt + 2 * y * radius  * cos( M_PI/6.0);
    auto c = start_pt - 2 * y * radius  * cos( M_PI/6.0);
    auto b = center -  2 * x * radius ;

    auto p = lerp(0.5,a,c);
    auto q = lerp(0.5,a,b);
    auto r = lerp(0.5,c,b);

    double weights[] = {1,0.5,1,0.5,1,0.5,1};
    Point cpts[] = {p, a   ,q    ,b    ,r   ,c, p};
    double ts[]  = {0, 0, 0, 2*M_PI/3,2*M_PI/3,4*M_PI/3,4*M_PI/3, 2*M_PI,2*M_PI,2*M_PI};
    auto spl = make_bspline(
        interleave( mk_stdvec(cpts,cpts+ sizeof(cpts)/sizeof(Point)),
                    std::vector<double>( weights , weights + sizeof(weights)/sizeof(double) ) ),
        std::vector<double>( ts , ts + sizeof(ts)/sizeof(double) ),
        2 );
    double st[] = {-2*M_PI/3,0,0};
    double es[] = {2*M_PI, 2*M_PI, 8.0*M_PI/3 };
    rebase_at_start(spl, st ).swap(spl);
    rebase_at_end(spl,es).swap(spl);
    return make_rbspline(std::move(spl));
}



#include "point.hpp"
using geom::point2d_t;
using geom::point3d_t;
/*
  Local Variables:
  eval:(load-file "./scripts/temp.el")
  eval: (setq methods (list "make_rbspline_from_circle" "make_circle"
  "foot_param" ) )
  eval: (setq pt-types  (list  "point2d_t" "point3d_t"  ) )
  eval:(instantiate-templates "circle" "geom" (list ) (product methods
  pt-types ) )
  End:
*/

#include "circle_inst.inl"


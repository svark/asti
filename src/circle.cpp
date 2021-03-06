//-*- mode:c++ -*-
#include "circle.hpp"
#include "geom_exception.hpp"
#include "rational_bspline_cons.hpp"
#include "smat.hpp"

// use geom::method convention to enable the elisp script below to
// instantiate

using geom::circle;
using geom::regular_tag;
using geom::rational_bspline;

template <class Point>
double
geom::foot_param(const circle<Point> &c,
                 const Point& p)
{
    auto x = normalize(c.start() - c.center());
    auto y = normalize(c.ydir());

    auto v = (p - c.center());

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

    if(tol::small(d))
        throw geom_exception(circle_too_large);

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

//  ../src/media/circle2.png
template <class Point>
rational_bspline< Point, regular_tag >
geom::make_rbspline_from_circle(const circle<Point>& circ)
{
    auto start_pt  =  circ.start();
    auto center    =  circ.center();
	auto const & x = normalize(start_pt - center);
    typedef decltype(circ.plane_normal()) normal_t;
    decltype(x)  y = normalize(cross(circ.plane_normal(), normal_t(x) ));

    double radius =  circ.radius();

    auto a = start_pt;
	auto b = start_pt + y * radius;
    auto c = center - x * radius +  y * radius;
    auto d = center - x * radius ;
	auto e = center - x * radius - y * radius;
	auto f = center + x * radius - y * radius;

    double weights[] = {1,0.5,0.5,1,0.5,0.5,1};
    Point cpts[] = {a,b,c,d,e,f,a};
    double ts[]  = {0, 0, 0,1.0/4,1.0/2,1.0/2,3.0/4,1,1,1};
    auto spl = make_bspline(
        interleave( mk_stdvec(cpts, cpts + sizeof(cpts)/sizeof(Point)),
                    mk_stdvec(weights , weights + sizeof(weights)/sizeof(double) ) ),
        std::vector<double>( ts , ts + sizeof(ts)/sizeof(double) ),
        2 );
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

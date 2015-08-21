#include "point.hpp"
#include "conic.hpp"
#include "geom_exception.hpp"

namespace geom {

template <class Point>
conic_arc < Point >
make_conic_arc_non_parallel(Point p[3],
                            decltype(Point() -  Point()) v[2])
{
    auto cps = closest_points(make_line(p[0], v[0]),
                              make_line(p[2], v[1]));

    assert(tol::small(len(cps.first - cps.second)));

    auto s = p[1]; // shoulder point
    auto cpsq = closest_points(make_line_joining(p[0], p[2]),
                               make_line_joining(cps.first, s));

    auto const & q = cpsq.first;

    assert(tol::small(len(q - cpsq.second)));

    double a = len(q - p[0])/len(q - p[2]);

    double u = a/(1+a);

    auto p1 = cps.first;

    double w = ( (1 - u) * (1 - u) *
                 dot(s - p[0], p1 - s)
                 + u * u * dot(s - p[2], p1 - s))
        /(2*u*(1 - u) * sqlen(p1 - p[1]));

    Point cpts[] = {p[0], p1, p[2]};
    double weights[] = {1, w, 1};
    return conic_arc<Point>(cpts, weights);
}

template <class Point>
conic_arc < Point >
make_conic_arc_parallel(Point p[3],
                        const decltype(Point() - Point())& v)
{
    auto s = p[1]; // shoulder point

    auto cpsq = closest_points(make_line_joining(p[0], p[2]),
                               make_line(s,v));

    auto const &  q = cpsq.first;

    if(tol::not_small(len(q - cpsq.second)))
		throw geom_exception(degenerate_or_small_conic);

    assert(tol::small(len(cpsq.first - cpsq.second)));

    double  a  = len(q - p[0])/len(q - p[2]);
    double  u  = a/(1 + a);
    double  b  = ((1-u)*(1-u)+u * u)/(2 * u * (1 - u));
    auto    p1 = make_pt((s-q)*b);
    double  w  = ((1 - u) * (1 - u) *
                  dot(s - p[0],
                      p1 - s) + u * u *
                  dot(s - p[2],
                      p1 - s)) / (2*u*(1 - u) *
                                  sqlen(p1 - p[1]));
    Point  cpts[]    = {p[0], p1, p[2]};
    double weights[] = {1, w, 1};
    return conic_arc<Point>(cpts, weights);
}

geom_error_code_t
conic_arc_preconditions(point3d_t p[], vector3d_t v[])
{
    if(tol::small(len(v[0])) || tol::small(len(v[1])))
		return tangent_vectors_too_small;

    auto nrml = cross(p[1] - p[0], p[2] - p[1]);

    if(tol::small(len(nrml)))
        return degenerate_or_small_conic;

    if(tol::not_small(dot(nrml, v[0])))
        return (vectors_not_in_plane_of_points);

    if(tol::not_small(dot(nrml, v[1])))
        return (vectors_not_in_plane_of_points);

    return no_error;
};


conic_arc<point3d_t>
make_conic_arc(point3d_t p[3], vector3d_t v[2])
{
    typedef conic_arc< point3d_t >::pointw_t pointw_t;

    auto err = conic_arc_preconditions(p,v);
    if(err)
        throw geom_error_code_t(err);

    if(tol::eq(fabs(dot(v[0], v[1])),
                len(v[0]) * len(v[1])))
    { // tangents are parallel
        return make_conic_arc_parallel(p, v[0]);
    }
    return  make_conic_arc_non_parallel(p,v);
}

conic_arc<point2d_t>
make_conic_arc(point2d_t p[3], vector2d_t v[2])
{
    typedef conic_arc< point3d_t >::pointw_t pointw_t;

    if( tol::eq(fabs(dot(v[0], v[1])),
                len(v[0]) * len(v[1])))
    { // tangents are parallel
        return  make_conic_arc_parallel(p, v[0]);
    }else
        return  make_conic_arc_non_parallel(p,v);
}


}

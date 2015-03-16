#include "conic.hpp"

namespace geom {

template <class Point>
conic_arc < Point >
make_conic_arc_non_parallel(Point p[3],
                            decltype(Point() - Point()) v[2])
{
    auto cps = closest_points(make_line(p[0], p[0] + v[0]),
                              make_line(p[2], p[2] + v[1]));

    if(!tol::eq(len(cps.first - cps.second), 0))
        throw geom_exception(lines_do_not_meet);

    auto s = p[1]; // shoulder point
    auto cpsq = closest_points(make_line(p[0], p[2]),
                               make_line(cps.first, s));

    auto const & q = cpsq.first;

    assert(tol::eq(len(q - cpsq.second),0));

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
    auto cpsq = closest_points(make_line(p[0], p[2]),
                               make_line(s,s+v));
    auto const &  q = cpsq.first;

    if(!tol::eq(len(q - cpsq.second), 0))
        throw geom_exception(lines_do_not_meet);

    assert(tol::eq(len(cpsq.first - cpsq.second),0));

    double a = len(q - p[0])/len(q - p[2]);

    double u = a/(1 + a);

    double b = ((1-u)*(1-u)+u * u)/(2 * u * (1 - u));

    auto p1 = make_pt((s-q)*b);

    double w = ( (1 - u) * (1 - u) * dot(s - p[0], p1 - s)
                 + u * u * dot(s - p[2], p1 - s))
        /(2*u*(1 - u) * sqlen(p1 - p[1]));

    Point cpts[] = {p[0], p1, p[2]};

    double weights[] = {1, w, 1};

    return conic_arc<Point>(cpts, weights);
}



conic_arc<point3d_t>
make_conic_arc(point3d_t p[3], vector3d_t v[2])
{
    typedef conic_arc< point3d_t >::pointw_t pointw_t;

    auto nrml = cross(p[1] - p[0], p[2] - p[1]);

    if(!tol::eq(dot(nrml, v[0]), 0.0))
        throw geom_exception(vectors_not_in_plane_of_points);
    if(!tol::eq(dot(nrml, v[1]), 0.0))
        throw geom_exception(vectors_not_in_plane_of_points);
    
    if( tol::eq(fabs(dot(v[0], v[1])),
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

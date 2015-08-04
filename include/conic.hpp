#ifndef ASTI_CONIC_HPP
#define ASTI_CONIC_HPP
#define _USE_MATH_DEFINES
#include <math.h>
#include "line.hpp"
#include "geom_exception.hpp"
#include "type_traits"

#include "tol.hpp"
#include "rational_bspline_cons.hpp"

namespace geom {

enum conic_type_t {parabola, hyperbola, ellipse} ;
//conic arc passing through three points start,intermediate point and
//an end
template <class Point>
struct conic_arc {
    typedef Point point_t;
    enum{ dim = point_dim <Point >::dimension };
    typedef typename inc_dimension<Point>::type pointw_t;

    conic_arc(point_t p_[],
              double weights[])
    {
        for(int i = 0;i < 3; ++i)
            if(tol::eq(weights[i], 0.0))
                p[i] = pointw_t(p_[i], 0.0);
            else {
                p[i] = pointw_t(p_[i], 1.0);
                p[i]*= weights[i];
            }
    }
    conic_arc(const pointw_t & p1,
              const pointw_t & p2,
              const pointw_t & p3)
    {
        p[0] = p1;
        p[1] = p2;
        p[2] = p3;
    }
    conic_arc(pointw_t p_[])
    {
        for(int i = 0;i < 3; ++i)
            p[i] = pointw_t(p_[i]);
    }

    double weight(int i) const { return p[i][dim]; }

    conic_type_t type() const {
        double w = weight(1);
        if(tol::eq(w,1))
            return parabola;
        else if(fabs(w) < 1)
            return ellipse;
        return hyperbola;
    }
    point_t eval(double u) const {
        pointw_t ep
            = make_pt((1 - u) * (1 - u) * make_vec( p[0] )
                      +  2 * u  * (1 - u) * make_vec( p[1] )
                      +  u * u  * make_vec( p[2] ) );
        point_t pep(ep);
        pep *= 1 / ep[dim];
        return pep;
    }

    std::tuple <pointw_t, pointw_t, pointw_t>
    split_conic_at_shoulder() const
    {
        // returns Q1, S, R1
        auto w = p[1][dim];
        auto q1 = lerp(0.5, p[0], p[1]);
        auto r1 = lerp(0.5, p[2], p[1]);
        auto s = lerp(0.5, q1, r1);
		s *= 1/s[dim];
		q1 *= sqrt((1 + w) / 2)/q1[dim]; r1 *= sqrt((1 + w) / 2)/r1[dim];
        return std::make_tuple( q1, s, r1);
    }
    pointw_t p[3];
};

template <class Point>
conic_arc < Point >
reverse_curve(const conic_arc < Point >& arc)
{
    return conic_arc<Point>(arc.p[2], arc.p[1], arc.p[0]);
}

template <class Point>
rational_bspline< Point, regular_tag >
make_rbspline_from_conic(const conic_arc<Point> &arc) {

    auto w = arc.weight(1);
    typedef typename inc_dimension<Point>::type PointW;
    enum{ dim = point_dim < Point >::dimension };
    if(w < -1.0)
    {
        auto car = reverse_curve(arc);
        car.p[1][dim] *=- 1;
        return make_rbspline_from_conic(car);
    }
    double alpha = angle_between(arc.p[1] - arc.p[0],
                                 arc.p[2] - arc.p[1]);
    if(w >= 1.0 || w > 0 && alpha > M_PI / 3) {
        // single segment parabola or hyperbola
        double ks[] =  {0, 0, 0, 1, 1, 1};
        auto spl(make_bspline_arr( arc.p, arc.p + 3, ks,
                                   ks+sizeof(ks)/sizeof(double), 2) );
        return make_rbspline(std::move(spl));
    }else if(w < 0 && alpha > M_PI / 2)
    {
        double ks[] = {0, 0, 0,
                       0.25, 0.25,
                       0.5, 0.5, 0.75, 0.75,
                       1, 1, 1};
        PointW q[3], s[3], r[3];
        std::tie(q[0], s[0], r[0]) = arc.split_conic_at_shoulder();
        conic_arc < Point > c1arc(arc.p[0], q[0], s[0]);
        conic_arc < Point > c2arc(s[0], r[0], arc.p[2]);
        std::tie(q[1], s[1], r[1]) = c1arc.split_conic_at_shoulder();
        std::tie(q[2], s[2], r[2]) = c2arc.split_conic_at_shoulder();
        PointW cpts[] = { arc.p[0],
                          q[1], s[1], r[1],
                          s[0],
                          q[2], s[2], r[2],
                          arc.p[2] };
        auto spl = make_bspline_arr(cpts,
                                    cpts + sizeof(cpts)/sizeof(PointW), ks,
                                    ks+sizeof(ks)/sizeof(double), 2);
        return make_rbspline(std::move(spl));
    }else
    {
        PointW q, s, r;
        std::tie(q, s, r) = arc.split_conic_at_shoulder();
        double ks[] = {0, 0, 0, 0.5, 0.5,1, 1, 1};
        PointW cpts[] = {arc.p[0],
                         q,s,r,
                         arc.p[2]};
        auto spl = make_bspline_arr(cpts, cpts+ sizeof(cpts)/sizeof(PointW),
                                    ks,ks+sizeof(ks)/sizeof(double), 2);
        return make_rbspline(std::move(spl));
    }
}

extern
conic_arc<point3d_t>
make_conic_arc(point3d_t p[3], vector3d_t v[2]);

extern
conic_arc<point2d_t>
make_conic_arc(point2d_t p[3], vector2d_t v[2]);


}
#endif // ASTI_CONIC_HPP

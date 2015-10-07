//-*- mode:c++ -*-
#include "line.hpp"
#include "point.hpp"

namespace geom {

point2d_t
intersect_lines(const line < point2d_t >& l1,
                const line < point2d_t >& l2)
{
    auto d1 = l1.direction();
    auto d2 = l2.direction();
    auto r =  l1.start_pt() - l2.start_pt();

    double b = dot(d1, d2);
    double f = dot(d2, r);

    double d = 1.0 - b * b;

    if(tol::eq(d, 0))
        return point2d_t::infinity();

    auto p1 = l1.start_pt() +  d1 * (dot(d1, (d2 * f - r )) / d);
    return p1;
}

}

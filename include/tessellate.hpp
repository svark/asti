#ifndef ASTI_TESSELLATE_HPP
#define ASTI_TESSELLATE_HPP

#include "special.hpp"
#include "polyline.hpp"
#include "spline_traits.hpp"
#include "bspline_queries.hpp"
#include "spline_norms.hpp"

namespace geom
{

template <class SplineCurve>
polyline<typename spline_traits<SplineCurve>::point_t>
tessellate(const SplineCurve& crv, double epsilon)
{
    double width = qry::param_range(crv).second;
    width -= qry::param_range(crv).first;
    assert(width > 0 );

	double twonrm =  two_norm_squared(crv);

    double delta =
        util::nroot( 24.0 * epsilon * epsilon * width
                     / twonrm, 5 );

    long num_segs = long (ceil(width/delta) );
    delta = width/num_segs;

    typedef typename spline_traits<SplineCurve>::point_t Point;
    decltype(mk_stdvec(Point())) crvpts(num_segs+1);
#pragma loop(hint_parallel(8))
    for(long i = 0 ; i <= num_segs;++i)
    {
        double u = crv.param_range().first + i*delta;
        crvpts[i] =  crv.eval(u);
    }
    return polyline<Point>(std::move(crvpts),
                           crv.param_range().first,
                           crv.param_range().second);
}
}
#endif // ASTI_TESSELLATE_HPP

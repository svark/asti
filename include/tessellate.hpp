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
	using geom::qry::start_param;
	using geom::qry::end_param;
	double s     = start_param(crv);
    double width = end_param(crv) - s;
    assert(width > 0 );
	
    double twonrmsq =  two_norm_squared(2,crv);

	double del = 24.0 *  width / twonrmsq;
    double delta =  sqrt(sqrt(del)) * sqrt(epsilon);
    long num_segs = long (ceil(width/delta) );
    delta = width/num_segs;

    typedef typename spline_traits<SplineCurve>::point_t Point;
    ARRAY_TYPE(Point) crvpts(num_segs+1);
#pragma loop(hint_parallel(8))
    for(long i = 0 ; i <= num_segs;++i)
    {
        double u = s + i*delta;
        crvpts[i] =  crv.eval(u);
    }
	
    return polyline<Point>(std::move(crvpts),
                           start_param(crv),
                           end_param(crv));
}
}
#endif // ASTI_TESSELLATE_HPP

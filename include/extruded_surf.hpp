#ifndef ASTI_EXTRUDED_SURF
#define ASTI_EXTRUDED_SURF

#include "bspline.hpp"
#include "bspline_surface.hpp"
#include "bspline_queries.hpp"
#include "type_utils.hpp"

namespace geom {

template <class SplineType>
bspline_surface<point3d_t,
                bspline_surface_traits<
                    typename spline_traits<SplineType>::rtag,
                    typename spline_traits<SplineType>::ptag,
                    geom::regular_tag,
                    extruded_surf >
                >
make_bspline_surface(
    const SplineType& curvU,
    const vector3d_t& dir
    )
{
	using namespace ops;
    typedef typename spline_traits<SplineType>::point_t point_t;
	typedef typename std::conditional<is_rational_type<SplineType>::value,point4d_t, point3d_t>::type pointw_t;
	typedef typename std::conditional<is_rational_type<SplineType>::value,vector4d_t, vector3d_t>::type vectorw_t;

	typedef bspline_surface<point3d_t,  bspline_surface_traits<
                    typename spline_traits<SplineType>::rtag,
                    typename spline_traits<SplineType>::ptag,
                    geom::regular_tag,
                    extruded_surf >> extruded_surf_t;

    typedef typename extruded_surf_t::cpts_t cpts_t;
    size_t stride = curvU.control_points().size();
    cpts_t pts(stride*2);
    
	static const std::integral_constant<int,point_t::dimension> spaceDim;

	typedef typename spline_traits<SplineType>::rtag rtag;
	for( auto pt : curvU.control_points() )
	{
		pointw_t p2;
		lift_dim(pt,rtag(),p2);
		pts.push_back(p2);
	}
    for( auto pt : curvU.control_points() )
    {
		pointw_t p2;
		vectorw_t dir2(dir,1.0);
		lift_dim(pt,rtag(),p2);
        pts.push_back(p2+dir2);
    }
    const double ts[] = {0,0,1,1};
    std::vector<double> ksv(ts,ts+sizeof(ts)/sizeof(double));
    return extruded_surf_t (
        pts,
        stride,
        curvU.knots(),
        ksv,
        curvU.degree(),
        1 );
}


}

#endif // ASTI_EXTRUDED_SURF

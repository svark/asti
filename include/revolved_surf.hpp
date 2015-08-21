#ifndef ASTI_REVOLVED_SURF
#define ASTI_REVOLVED_SURF

#include "point_fwd.hpp"
#include "line.hpp"
#include "spline_traits.hpp"
#include "bspline_fwd.hpp"
#include "bspline_surface.hpp"
#include "conic.hpp"
namespace geom{

template <class T> struct line;
template <class SplineCurve>
auto make_revolved_bspline_surf(
    const SplineCurve & spl,
    const line<point3d_t> &axis,
    double theta) -> bspline_surface<point3d_t,
                                     bspline_surface_traits<
                                         rational_tag,
                                         typename spline_traits<SplineCurve>::ptag,
                                         regular_tag, revolved_surf> > 
{
    typedef bspline_surface<point3d_t,
                            bspline_surface_traits<
                                rational_tag,
                                typename spline_traits<SplineCurve>::ptag,
                                regular_tag, revolved_surf> >
        surf_t;
    auto const &cpts = spl.control_points();
    auto  stride = cpts.size();
    double ct2 = cos(theta/2);
    double st2 = sin(theta/2);
    double ct = cos(theta);
    double st = sin(theta);
	typedef decltype(mk_stdvec(point4d_t())) arr_t;
    std::unique_ptr<arr_t> pts;
    std::vector<double> t_v;
    typename spline_traits<SplineCurve>::rtag rtagu;

	using namespace qry;
	size_t i = 0;
    for(auto pt : cpts)
    {
        double w    = weight(pt,rtagu);
        bool isdir  = tol::eq(w,0);
        auto p3d    = point3d_t(pt);
        auto onaxis = point3d_t(0.0);
        if(!isdir) { //not points at infinity
            p3d = auto_lift_dim3(pt,rtagu,polynomial_tag());
			onaxis = closest_point_on_line(axis,p3d);
        }

        auto radial_dir     = p3d - onaxis;
		auto ydir           = cross(radial_dir,axis.direction());
		auto tangential_dir = tol::not_small(len(ydir))? normalize(ydir) : vector3d_t(0.0);
        auto mid            = onaxis + ct2*radial_dir + st2*tangential_dir;
        auto end            = onaxis + ct*radial_dir  + st*tangential_dir;

        point3d_t vs[] = {p3d,mid,end};

        auto const & arc =
            make_rbspline_from_conic(make_circular_arc(vs));
        if(!t_v.size())
            t_v = arc.knots();
        int ncpts = num_cpts(arc);
        if(!pts)
            pts.reset(new arr_t(ncpts*stride));

        for(int j = 0; j < ncpts; ++j)
        {
            double wm = weight(arc.control_points()[j], rational_tag());
            (*pts)[i+j*stride] = point4d_t(point3d_t(arc.control_points()[j]), wm*w);
        }
		++i;
    }

    return surf_t(*pts,
                  cpts.size(),
                  spl.knots(),
                  t_v,
				  spl.degree(), 2
        );
}
/*
template <class SplineCurve>
bspline_surface<typename SplineCurve::point_t,
                bspline_surface_traits<
                    rational_tag,
                    typename spline_traits<SplineCurve>::ptag,
                    periodic_tag, revolved_surf> >
make_revolved_bspline_surf(
    const SplineCurve & spl,
    const line<point3d_t> &axis
    )
{
    static_assert(point_dim<typename SplineCurve::point_t>::dimension==3,
                  "expecting a 3 dimensional spline curve");

    typedef typename SplineCurve::point_t pt_t;
    typedef bspline_surface<pt_t,
                            bspline_surface_traits<
                                rational_tag,
                                typename spline_traits<SplineCurve>::ptag,
                                periodic_tag, revolved_surf> >
        surf_t;

    typedef typename spline_traits<SplineCurve>::rtag rtag;

    auto const &cpts = spl.control_points();
    std::vector<double> t_v;

    std::unique_ptr<std::vector<point4d_t>> pts;
    int i = 0;
    for(auto pt : cpts)
    {
        auto const &center = closest_point_on_line(line, point3d_t(pt));
        auto const & splv  =
            make_rbspline_from_circle(
                circle<decltype(center)>( center, point3d_t(pt),
                                          line.direction() ));
        if(!t_v.size())
            t_v = splv.knots();

        int ncpts = num_cpts(splv);

        if(!pts)
            pts.reset(new std::vector<point4d_t>[ncpts*stride]);

        int j = 0;
        for(auto cptv : splv.control_points())
        {
            double wu = weight(pt,rtag, std::integral_constant<int,3>());
            pts[i+j*stride] = cptv;
            pts[i+j*stride]*= wu;
            ++j;
        }
        ++i;
    }

    return surf_t(std::move(*pts),
                  cpts.size(),
                  spl.knots(),
                  t_v
        );
        } */
}

#endif // ASTI_REVOLVED_SURF

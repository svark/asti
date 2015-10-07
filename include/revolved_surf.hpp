#ifndef ASTI_REVOLVED_SURF
#define ASTI_REVOLVED_SURF

#include "point_fwd.hpp"
#include "line.hpp"
#include "spline_traits.hpp"
#include "bspline_fwd.hpp"
#include "bspline_surface.hpp"
#include "conic.hpp"
#include "circle.hpp"
namespace geom{

namespace impl {
template <class SplineCurve,class PTag, class MakeFn>
auto make_revolved_bspline_surf(
    const SplineCurve & spl,
    const line<point3d_t> &axis,
    double theta, PTag tag, MakeFn mk_conic_fn )
    -> bspline_surface<point3d_t,
    bspline_surface_traits<
    rational_tag,
    typename spline_traits<SplineCurve>::ptag,
    PTag, revolved_surf> >

{
    typedef bspline_surface<point3d_t,
                            bspline_surface_traits<
                                rational_tag,
                                typename spline_traits<SplineCurve>::ptag,
                                PTag, revolved_surf> >
        surf_t;
    auto const &cpts = spl.control_points();
    auto  stride = cpts.size();
    assert(stride != 0);
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
        if(!isdir) { //not a point at infinity
            onaxis = closest_point_on_line(axis,p3d);
        }

        auto radial_dir     = p3d - onaxis;
        auto ydir           = cross(radial_dir,axis.direction());
        bool ydir_notsmall  = tol::not_small(len(ydir));
        auto tangential_dir = ydir_notsmall? normalize(ydir)  : vector3d_t(0.0);
        auto mid            = onaxis + ct2*radial_dir + st2*tangential_dir;
        auto end            = onaxis + ct*radial_dir  + st*tangential_dir;

        const point3d_t vs[] = {p3d, mid, end};

        auto const & arc = mk_conic_fn(vs);

        if(!t_v.size())
            t_v = arc.knots();

        int ncpts = num_cpts(arc);
        if(!pts)
            pts.reset(new arr_t(ncpts*stride));

        for(int j = 0; j < ncpts; ++j)
        {
            double wm = weight(arc.control_points()[j], rational_tag());
            (*pts)[i+j*stride] = point4d_t(point3d_t(
                                               arc.control_points()[j]),
                                           wm*w);
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

}
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
    auto make_fn = [] (const point3d_t* ps) {
        return make_rbspline_from_conic(make_circular_arc(ps));
    };
    return impl::make_revolved_bspline_surf(spl,axis,theta,
                                            regular_tag(),
                                            make_fn);
}

template <class SplineCurve>
auto make_revolved_bspline_surf(
    const SplineCurve & spl,
    const line<point3d_t> &axis)
    -> bspline_surface<point3d_t,
    bspline_surface_traits<
    rational_tag,
    typename spline_traits<SplineCurve>::ptag,
    periodic_tag, revolved_surf> >
{
    auto make_fn = [] (const point3d_t* ps) {
        return make_rbspline_from_circle(make_circle(ps[0],ps[1],ps[2]));
    };
    return impl::make_revolved_bspline_surf(spl,
                                            axis,M_PI/*this is rather
                                                      * arbitrary*/,
                                            periodic_tag(),
                                            make_fn);
}

}
#endif // ASTI_REVOLVED_SURF

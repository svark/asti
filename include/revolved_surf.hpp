#ifndef ASTI_REVOLVED_SURF
#define ASTI_REVOLVED_SURF
#include "point.hpp"
#include "line.hpp"
#if 0
namespace geom{

template <class SplineCurve>
bspline_surface<typename SplineCurve::point_t,
                bspline_surface_traits<
                    rational_tag,
                    typename spline_traits<SplineCurve>::ptag,
                    regular_tag, revolved_surf> >
make_revolved_bspline_surf(
    const SplineCurve & spl,
    const line<point3d_t> &axis,
    double theta)
{
    typedef typename SplineCurve::point_t pt_t;
    typedef bspline_surface<pt_t,
                            bspline_surface_traits<
                                rational_tag,
                                typename spline_traits<SplineCurve>::ptag,
                                regular_tag, revolved_surf> >
        surf_t;
    auto const &cpts = spl.control_points();
    double ct2 = cos(theta/2);
    double st2 = sin(theta/2);
    double ct = cos(theta);
    double st = sin(theta);
    std::unique_ptr<std::vector<point4d_t>> pts;
    std::vector<double> t_v;

    for(auto pt : cpts)
    {
        auto onaxis = closest_point_on_line(line, pt);
        auto radial_dir     = pt - onaxis;
        auto tangential_dir = normalize(cross(radial_dir,axis));
        auto m = onaxis + ct2*radial_dir + st2*tangential_dir;
        auto q = onaxis + ct*radial_dir  + st*tangential_dir;

        decltype(m) vs[] = {pt,m,q};

        auto const & spl =
            make_rbspline_from_conic(make_circular_arc(vs));
        if(!t_v.size())
            t_v = spl.knots();
        if(!pts)
            pts.reset(new std::vector<point4d_t>[ncpts*stride]);
        int ncpts = num_cpts(spl);

        for(int j = 0; j < numcpts; ++j)
        {
            pts[i+j*stride] = spl.control_points()[j];
        }
    }

    return surf_t(flatten_cpts(pts),
                  cpts.size(),
                  spl.knots(),
                  t_v
        );
}

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
}
#endif

#endif // ASTI_REVOLVED_SURF

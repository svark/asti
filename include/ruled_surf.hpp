#ifndef ASTI_RULED_SURF
#define ASTI_RULED_SURF
#include "raise_degree.hpp"
#include "insert_knot.hpp"
#include "bspline_surface.hpp"
#include "bspline_queries.hpp"

namespace geom
{

template <class SplineType>
bspline_surface < typename SplineType::point_t,
                  bspline_surface_traits <
                      typename spline_traits<SplineType>::rtag,
                      typename spline_traits<SplineType>::ptag,
                      regular_tag, ruled_surf >  >
make_ruled_bspline_surface(
    const SplineType & spl1,
    const SplineType & spl2)
{
    auto s1(spl1);
    auto s2(spl2);

    using namespace ops;
    match_degrees(s1, s2);
    match_knots(s1, s2);

    typedef typename SplineType::point_t point_t;
    typename SplineType::cpts_t surf_cpts;

    size_t stride = num_cpts(s1);
    assert(stride == num_cpts(s2));

    surf_cpts.reserve(2 * stride);
    surf_cpts.assign(s1.control_points().begin(), s1.control_points().end());
    surf_cpts.insert(surf_cpts.end(),
                     s2.control_points().begin(),
                     s2.control_points().end());
    double ks[] = {0, 0, 1, 1};
    std::vector<double> ksv(ks, ks + sizeof(ks) / sizeof(double));

    typedef  bspline_surface_traits < typename spline_traits<SplineType>::rtag,
                                      typename spline_traits<SplineType>::ptag,
                                      regular_tag, ruled_surf > traits_t;

    return bspline_surface<point_t,traits_t>
        (
            std::move(surf_cpts),
            stride,
            s1.knots(),
            ksv,
            s1.degree(),
            1
         );
}

}
#endif // ASTI_RULED_SURF

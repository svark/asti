#ifndef ASTI_RULED_SURF
#define ASTI_RULED_SURF
#include "smat.hpp"
#include "raise_degree.hpp"
#include "insert_knot.hpp"
namespace geom
{

template <class SplineType>
struct ruled_bspline_surface {

    typedef typename SplineType::point_t point_t;
    typedef typename SplineType::vector_t vector_t;

    ruled_bspline_surface(const SplineType & spl1_,
                          const SplineType & spl2_)
        :spl1(spl1_), spl2(spl2_),
         surf(make_ruled_bspline_surface(spl1, spl2))
    {
    }
    point_t eval(double u,  double v)
    {
        return surf.eval(u, v);
    }

    SplineType spl1, spl2;
    bspline_surface < point_t > surf;
};


template <class SplineType>
bspline_surface < point_t,
                  bspline_surface_traits <
                      typename spline_traits<SplineType>::rtag,
                      typename spline_traits<SplineType>::ptag,
                      regular_tag >  >
make_ruled_bspline_surface(
    const SplineType & spl1,
    const SplineType & spl2)
{
    auto s1(extract_regular_curve(spl1).deoptimize());
    auto s2(extract_regular_curve(spl2).deoptimize());

    match_degrees(s1, s2);
    match_knots(s1, s2);

    typedef typename SplineType::point_t point_t;
    std::vector<point_t> surf_cpts;
    uint stride =  s1.control_points().size();
    surf_cpts.reserve(2 * stride);
    surf_cpts.assign(s1.control_points());
    surf_cpts.insert(surf_cpts.end(), s2.control_points());
    double ks[] = {0, 0, 1, 1};
    std::vector<double> ksv(ks, ks + sizeof(ks) / sizeof(double));
    return bspline_surface < point_t,
                             bspline_surface_traits <
                                 typename spline_traits<SplineType>::rtag,
                                 typename spline_traits<SplineType>::ptag,
                                 regular_tag >  >
        (
            std::move(surf_cpts),
            stride,
            s.knots(),
            ksv
            );
}


}
#endif // ASTI_RULED_SURF

#ifndef BSPLINE_X_CONS
#define BSPLINE_X_CONS
#include "spline_traits.hpp"
#include "bspline_cons.hpp"
#include "rational_bspline_cons.hpp"
#include "periodic_bspline_cons.hpp"

namespace geom {
namespace impl {
template<class SplineType, class CptsT, class KnotsT>
SplineType make_bsplinex(CptsT && pts, KnotsT && ks, int deg, regular_tag, nonrational_tag)
{
    return geom::make_bspline(std::forward < CptsT > (pts),
                              std::forward < KnotsT > (ks), deg);
}

template<class SplineType, class CptsT, class KnotsT>
SplineType make_bsplinex(CptsT && pts, KnotsT && ks, int deg,
                         periodic_tag, nonrational_tag)
{
    return geom::make_periodic_bspline(std::forward < CptsT > (pts),
                                       std::forward < KnotsT > (ks),
                                       deg);
}

template<class SplineType, class CptsT, class KnotsT>
SplineType make_bsplinex(CptsT && pts, KnotsT && ks, int deg,
                         regular_tag, rational_tag)
{
    return make_rbspline(
        geom::make_bspline(std::forward < CptsT > (pts),
                           std::forward < KnotsT > (ks), deg));
}

template<class SplineType, class CptsT, class KnotsT>
SplineType make_bsplinex(CptsT && pts, KnotsT && ks, int deg,
                         periodic_tag, rational_tag)
{
    return make_rbspline(
        geom::make_periodic_bspline(std::forward < CptsT > (pts),
                                    std::forward < KnotsT > (ks), deg));
}
}

template<class SplineType, class CptsT, class KnotsT >
SplineType
make_bsplinex(CptsT pts, KnotsT ks, int deg)
{
    return impl::make_bsplinex < SplineType > (
        std::move(pts), std::move(ks),
        deg,
        typename spline_traits<SplineType >::ptag(),
        typename spline_traits<SplineType >::rtag());
}

}
#endif // BSPLINE_X_CONS

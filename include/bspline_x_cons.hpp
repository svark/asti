#ifndef ASTI_ASTI_BSPLINE_X_CONS
#define ASTI_ASTI_BSPLINE_X_CONS
#include "spline_traits.hpp"
#include "bspline_cons.hpp"
#include "rational_bspline_cons.hpp"
#include "periodic_bspline_cons.hpp"

namespace geom {

template<class CptsT, class KnotsT>
typename get_traits_type_from_tags < regular_tag,  polynomial_tag >::type::spline_type
make_bspline(CptsT && pts, KnotsT && ks, int deg,
             regular_tag, polynomial_tag)
{
    return geom::make_bspline(std::forward < CptsT > (pts),
                              std::forward < KnotsT > (ks), deg);
}

template<class CptsT, class KnotsT>
typename get_traits_type_from_tags < periodic_tag,  polynomial_tag >::type::spline_type
make_bspline(CptsT && pts, KnotsT && ks, int deg,
             periodic_tag, polynomial_tag)
{
    return geom::make_periodic_bspline(std::forward < CptsT > (pts),
                                       std::forward < KnotsT > (ks),
                                       deg);
}

template< class CptsT, class KnotsT>
typename get_traits_type_from_tags < regular_tag,  rational_tag >::type::spline_type
make_bspline(CptsT && pts, KnotsT && ks, int deg,
             regular_tag, rational_tag)
{
    return make_rbspline(
        geom::make_bspline(std::forward < CptsT > (pts),
                           std::forward < KnotsT > (ks), deg));
}

template<class CptsT, class KnotsT>
typename get_traits_type_from_tags < periodic_tag,  rational_tag >::type::spline_type
make_bspline(CptsT && pts, KnotsT && ks, int deg,
             periodic_tag, rational_tag)
{
    return make_rbspline(
        geom::make_periodic_bspline(std::forward < CptsT > (pts),
                                    std::forward < KnotsT > (ks), deg));
}

template<class CptsIterT,
         class KnotsIterT,
         class PTag,
         class RTag>
typename get_traits_type_from_tags <PTag, RTag>::type::spline_type
make_bspline(CptsIterT  ptsf,
             CptsIterT  ptsl,
             KnotsIterT ksf,
             KnotsIterT ksl,
             int deg,
             PTag, RTag)
{
    return make_bspline(mkstdvec(ptsf, ptsl),
                        mkstdvec(ksf, ksl),
                        deg, PTag(), RTag());
}

}
#endif // ASTI_BSPLINE_X_CONS

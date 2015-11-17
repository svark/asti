#ifndef ASTI_ASTI_BSPLINE_X_CONS
#define ASTI_ASTI_BSPLINE_X_CONS
#include "spline_traits.hpp"
#include "bspline_cons.hpp"
#include "rational_bspline_cons.hpp"
#include "periodic_bspline_cons.hpp"

namespace geom {

template<class CptsT, class KnotsT>
auto
make_bspline(CptsT && pts, KnotsT && ks, int deg,
             regular_tag, polynomial_tag)
{
    return make_bspline(std::forward < CptsT > (pts),
                        std::forward < KnotsT > (ks), deg);
}

template<class CptsT, class KnotsT>
auto
make_bspline(CptsT && pts, KnotsT && ks, int deg,
             periodic_tag, polynomial_tag)
{
    return make_periodic_bspline(std::forward<CptsT>(pts),
                                 std::forward<KnotsT> (ks),
                                 deg);
}

template<class CptsT, class KnotsT,class PTag>
auto
make_bspline(CptsT && pts, KnotsT && ks, int deg,
             PTag , rational_tag)
{
    return make_rbspline(std::move(make_bspline(pts,ks,deg,PTag(),
                                                polynomial_tag())));
}

}
#endif // ASTI_BSPLINE_X_CONS

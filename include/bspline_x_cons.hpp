#ifndef ASTI_ASTI_BSPLINE_X_CONS
#define ASTI_ASTI_BSPLINE_X_CONS
#include "spline_traits.hpp"
#include "bspline_cons.hpp"
#include "rational_bspline_cons.hpp"
#include "periodic_bspline_cons.hpp"

namespace geom {

// to remove when you have c++14
#define MAKE_BSPLINE  make_bspline(std::forward < CptsT > (pts),        \
                                   std::forward < KnotsT > (ks), deg)

#define MAKE_PERIODIC_BSPLINE   make_periodic_bspline(std::forward<CptsT> (pts), \
                                                      std::forward<KnotsT> (ks), \
                                                      deg)
#define MAKE_RBSPLINE make_rbspline(std::move(make_bspline(pts,ks,deg,PTag(),polynomial_tag())))

template<class CptsT, class KnotsT>
    auto
    make_bspline(CptsT && pts, KnotsT && ks, int deg,
                 regular_tag, polynomial_tag) ->
    RAWTYPE(MAKE_BSPLINE)
{
    return MAKE_BSPLINE;
}

template<class CptsT, class KnotsT>
auto
make_bspline(CptsT && pts, KnotsT && ks, int deg,
             periodic_tag, polynomial_tag) ->
    RAWTYPE(MAKE_PERIODIC_BSPLINE)
{
    return MAKE_PERIODIC_BSPLINE;
}

template<class CptsT, class KnotsT,class PTag>
auto
make_bspline(CptsT && pts, KnotsT && ks, int deg,
             PTag , rational_tag) ->
    RAWTYPE(MAKE_RBSPLINE)
{
    return MAKE_RBSPLINE;
}

/*template<class CptsIterT,
  class KnotsIterT,
  class PTag,
  class RTag>
  auto
  make_bspline(CptsIterT  ptsf,
  CptsIterT  ptsl,
  KnotsIterT ksf,
  KnotsIterT ksl,
  int deg,
  PTag, RTag) -> RAWTYPE( make_bspline(mkstdvec(ptsf, ptsl),
  mkstdvec(ksf, ksl),
  deg, PTag(), RTag()))
  {
  return make_bspline(mkstdvec(ptsf, ptsl),
  mkstdvec(ksf, ksl),
  deg, PTag(), RTag());
  }*/

}
#endif // ASTI_BSPLINE_X_CONS

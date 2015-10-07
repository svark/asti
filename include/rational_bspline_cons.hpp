#ifndef ASTI_RATIONAL_BSPLINE_CONS
#define ASTI_RATIONAL_BSPLINE_CONS

#include "periodic_bspline_cons.hpp"
#include "rational_bspline.hpp"
#include "tol.hpp"

namespace geom {

template<class PointContT, class WeightsT>
auto interleave(const PointContT & ps,
                const WeightsT & ws)
    -> RAWTYPE(mk_stdvec(higher_dim(ps[0])))
{
    typedef typename inc_dimension < RAWTYPE(ps[0]) >::type pointw_t;
    ARRAY_TYPE(pointw_t) pws;
    pws.reserve(ps.size());
    auto witer = std::begin(ws);
    enum {dim = point_dim<pointw_t>::dimension };
    for(const auto &p: ps)
    {
        pointw_t pw(p,1.0);
        if(tol::eq(*witer, 0))
            pw[dim] = 0.0;
        else
            pw *= *witer;
        ++witer;
        pws.push_back(pw);
    }
    return pws;
}

template <class CptsT, class WeightsT>
auto
make_periodic_rbspline(const CptsT & pts,
                       const WeightsT & weights,
                       std::vector<double> ks,
                       int degree_) -> rational_bspline<RAWTYPE(pts[0]), periodic_tag>
{
    return make_rbspline(make_periodic_bspline_wrap (
                             std::move(interleave(pts, weights) ),
                             std::move(ks), degree_));
}

template <class CptsT>
auto
make_periodic_rbspline(CptsT  pts,
                       std::vector<double> ks,
                       int degree_) -> rational_bspline<RAWTYPE(lower_dim(pts[0])), periodic_tag>
{
    return make_rbspline(make_periodic_bspline_wrap(
                             std::move(pts),
                             std::move(ks), degree_));
}

template <class CptsT, class WeightsT>
auto
make_rbspline(const CptsT & pts,
              const WeightsT & weights,
              std::vector<double> ks,
              int degree_) -> rational_bspline<RAWTYPE(pts[0]), regular_tag>
{
    return make_rbspline(make_bspline (
                             std::move(interleave(pts, weights) ),
                             std::move(ks), degree_));
}

}
#endif // ASTI_RATIONAL_BSPLINE_CONS

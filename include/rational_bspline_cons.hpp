#ifndef ASTI_RATIONAL_BSPLINE_CONS
#define ASTI_RATIONAL_BSPLINE_CONS

#include "periodic_bspline_cons.hpp"
#include "rational_bspline.hpp"

namespace geom {

template<int dim, class PointContT, class WeightsT>
std::vector<typename inc_dimension < pt_t<dim> >::type,
            typename point_dim<typename inc_dimension < pt_t<dim> >::type >::alloc_t>
interleave(const PointContT & ps,
           const WeightsT & ws)
{
    typedef typename inc_dimension < pt_t<dim> >::type pointw_t;
    std::vector<pointw_t, typename point_dim<pointw_t>::alloc_t > pws;
    pws.reserve(ps.size());
    auto witer = std::begin(ws);
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
rational_bspline<periodic_bspline <RAWTYPE(pts[0])>>
make_periodic_rbspline(const CptsT & pts,
                       const WeightsT & weights,
                       std::vector<double> ks,
                       int degree_)
{
    return make_rbspline(make_periodic_bspline ( spline_wrap_t(),
                         std::move(interleave<dimension>(pts, weights) ),
                         std::move(ks), degree_));
}

template <class CptsT, class WeightsT>
rational_bspline<bspline<RAWTYPE(pts[0])>>
make_rbspline(const CptsT & pts,
              const WeightsT & weights,
              std::vector<double> ks,
              int degree_)
{
    return make_rbspline(make_bspline (
                         std::move(interleave<dimension>(pts, weights) ),
                         std::move(ks), degree_));
}

}
#endif // ASTI_RATIONAL_BSPLINE_CONS

#ifndef ASTI_RATIONAL_BSPLINE_CONV
#define ASTI_RATIONAL_BSPLINE_CONV

#include "periodic_bspline_cons.hpp"
#include "rational_bspline.hpp"

namespace geom {

template <class CptsT, class WeightsT>
rational_bspline<RAWTYPE(pts[0]), periodic_tag>
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
rational_bspline<RAWTYPE(pts[0]), periodic_tag>
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
#endif // ASTI_RATIONAL_BSPLINE_CONV

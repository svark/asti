#ifndef ASTI_SPLINE_NORMS_HPP
#define ASTI_SPLINE_NORMS_HPP
#include "spline_traits.hpp"
#include "point.hpp"
#include "combinations.hpp"
#include "split_into_bezier_patches.hpp"
#include "hodograph.hpp"
#include "integrate1d.hpp"

namespace geom {
namespace impl {
template <class SplineType>
double two_norm_squared(const SplineType &spl, rational_tag)
{
    auto const &crvs    =
        ops::split_into_bezier_patches(
            qry::get_spline(spl));

    double _2norm = 0;
    for(auto c : crvs ) {
        auto accel = [&c] (double t)
            {
                return sqlen(c.eval_derivative(2,t));
            };
        _2norm += integral(accel,
                           c.param_range().first,
                           c.param_range().second );
    }
    return _2norm;
}

template <class SplineType>
double two_norm_squared(const SplineType &spl, polynomial_tag)
{
    int         n       = spl.degree() - 2;
    auto const  &ncis   = util::ncks(n);
    auto const  &_2ncis = util::ncks(2*n);
    auto const  &crvs   =
        ops::split_into_bezier_patches(
            qry::hodograph(qry::get_spline(spl),2));

    double _2norm = 0.0;
    for(auto c : crvs)
    {
        double thisnorm = 0;
        for (int i = 0; i <= n; ++i) {
            size_t ni = ncis[i];
            for (int j = 0; j <= n; ++j) {
                thisnorm +=  ni * ncis[j] *
                    dot( make_vec(c.control_points()[i]),
                         make_vec(c.control_points()[j]))/
                    _2ncis[i+j];
            }
        }
        double w = c.param_range().second;
        w-= c.param_range().first;
        _2norm += (thisnorm*w);
    }
    return _2norm/(2*n+1);
}
}

template <class SplineType>
double two_norm_squared(const SplineType &spl)
{
    return impl::two_norm_squared(spl, typename spline_traits<SplineType>::rtag());
}

}
#endif // ASTI_SPLINE_NORMS_HPP

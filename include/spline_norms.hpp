#ifndef ASTI_SPLINE_NORMS_HPP
#define ASTI_SPLINE_NORMS_HPP
#include "spline_traits.hpp"
#include "point.hpp"
#include "combinations.hpp"
#include "split_into_bezier_patches.hpp"
#include "hodograph.hpp"
#include "integrate1d.hpp"
#include "bspline_queries.hpp"

namespace geom {
namespace impl {
template <class SplineType>
double two_norm_squared(int der, const SplineType &spl, rational_tag)
{
    double _2norm = 0;
	using namespace geom::qry;

    auto accel = [&spl,der] (double t)
        {
			auto sder = derivative_approx(spl,t,der);
            return sqlen(sder);
        };
    _2norm += integral(accel,
                        start_param(spl),
                        end_param(spl));
    return _2norm;
}

template <class SplineType>
double two_norm_squared(int der, const SplineType &spl, polynomial_tag)
{
    int         n       = spl.degree() - 2;
    auto const  &ncis   = util::ncks(n);
    auto const  &_2ncis = util::ncks(2*n);
    auto const  &crvs   =
        ops::split_into_bezier_patches(
            qry::hodograph(qry::get_spline(spl),der));

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
		double s,e;
		std::tie(s,e) = c.param_range();
        double w = e - s;
        _2norm += (thisnorm*w);
    }
    return _2norm/(2*n+1);
}
}

template <class SplineType>
double two_norm_squared(int der, const SplineType &spl)
{
    return impl::two_norm_squared(der, spl, 
		typename spline_traits<SplineType>::rtag());
}

}
#endif // ASTI_SPLINE_NORMS_HPP

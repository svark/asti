#include "periodic_bspline_cons.hpp"
#include <numeric>
namespace geom {
namespace ops {

template <class Point>
periodic_bspline < Point >
rotate_base_knot(const periodic_bspline < Point >& pspl, size_t nu)
{
    int p = pspl.degree();
    auto const & cpts = pspl.control_points();
    auto const &t = pspl.knots();
    typename periodic_bspline < Point >::cpts_t
        cpts_unwrapped(cpts.cbegin() + p, cpts.cend());
    std::rotate(cpts_unwrapped.begin(),
                cpts_unwrapped.begin() + nu - p,
                cpts_unwrapped.end() );
    std::vector<double> knots_unwrapped;
    knots_unwrapped.reserve (t.size());
    std::rotate_copy(t.begin() + p, t.begin() + nu, t.end() - p,
                     std::back_inserter(knots_unwrapped));
    std::vector<double> startDiffs(nu - p + 1);
    std::adjacent_difference(t.begin() + p, t.begin() + nu + 1,
                             startDiffs.begin());
    startDiffs[0] =  t.end()[-p-1];
    std::partial_sum(startDiffs.begin(), startDiffs.end(),
                     knots_unwrapped.end() - nu + p - 1);
    return make_periodic_bspline_wrap(cpts_unwrapped, knots_unwrapped, p);
}

}
}

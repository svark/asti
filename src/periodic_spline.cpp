#include "stdafx.h"
#include "periodic_spline.hpp"
#include <algorithm>
#include <numeric>
namespace geom
{
template <class Point>
typename periodic_bspline<Point>::tuple_t
periodic_bspline<Point>::wrap( const cpts_t& pts,
                               const knots_t &ks,
                               int degree)
{
    assert(degree >= 1);
    assert(ks.size() && pts.size() && ks.size() == pts.size()+1);
    cpts_t cpts;
    cpts.reserve(pts.size() + degree);
    // within the range, t_\mu, t_{\mu+1} the spline evaluation
    // depends on cpt_\mu-d..\cpt_\mu
    cpts.assign(pts.end() - degree, pts.end());
    cpts.insert(cpts.end(),pts.begin(), pts.end());

    // within the range, t_\mu, t_{\mu+1} the spline evaluation
    // depends on t_{\mu-d+1}.....t_\mu,t_{\mu+1},....t_{\mu+d}
    // (check the R_d matrix) ./media/rmatrix.png at \mu = 0, we need to
    // ensure that t_{\mu-d+1} and accessible, so d-1 extra knots
    // are prefixed here at \mu = ks.size()-2, we need to ensure
    // that t_{\mu+d} is accessible, so d-1 extra knots are appended
    // here as well.
    knots_t t;
    t.reserve( ks.size() + 2*degree );
    knots_t startDiffs(degree), endDiffs(degree);

    // get the adjacent diffs between ks_0,ks_1,...
    std::adjacent_difference(ks.begin(), ks.begin() + degree,
                             startDiffs.begin());
    // get the adjacent diffs between ks_{n-1},ks_{n-2}, ....
    std::adjacent_difference(ks.rbegin(), ks.rbegin() + degree,
                             endDiffs.begin());


    startDiffs[0] = ks.back(); //change to ks_{n-1} in preparation
    //for the partial sums below
    endDiffs[0] = ks.front(); //change the start of backwardly
    //oriented end diffs to ks_0

    //  ./media/periodic.png

    std::vector<double> back_plus(startDiffs.size()),
        front_minus(endDiffs.size());

    std::partial_sum( startDiffs.begin(),
                      startDiffs.end(),
                      back_plus.begin() );

    std::partial_sum( endDiffs.begin(),
                      endDiffs.end(),
                      front_minus.begin() );

    t.assign(front_minus.rbegin(), front_minus.rend());
    t.insert(t.end(),ks.begin() + 1, ks.end()); // already copied
                                                // ks.front() so the +1
    t.insert(t.end(),back_plus.begin() + 1, back_plus.end());
                             //already copied ks.back()
    return std::forward_as_tuple(std::move(cpts), std::move(t),
                                 std::move(degree));
}

}
#include "point.hpp"
namespace geom {
template <class Point> struct periodic_bspline;
/*
Local Variables:
eval:(load-file "./temp.el")
eval:(instantiate-templates "periodic_bspline" '("double" "point2d_t" "point3d_t"  "point4d_t")  '() '() )
End:
*/
#include "periodic_spline_inst.cpp"
}
template struct geom::periodic_bspline<geom::pt_t<3>>;
template struct geom::periodic_bspline<geom::pt_t<4>>;
template struct geom::periodic_bspline<double>;

#ifndef ASTI_PERIODIC_BSPLINE_CONS
#define ASTI_PERIODIC_BSPLINE_CONS
#include "bspline_cons.hpp"
#include "periodic_bspline.hpp"
#include "type_utils.hpp"
#include <utility>
#include <numeric>
#include "bspline_queries.hpp"
namespace geom {

template <class CptsT, class KnotsT>
std::tuple<CptsT,std::vector<double>, int>
wrap_spline_params( const CptsT& pts,
                    const KnotsT& ks,
                    int degree)
{
    assert(degree >= 1);
    assert(pts.size() && ks.size() <= pts.size()+1);
    CptsT cpts;
    cpts.reserve(pts.size() + degree);
    // within the range, t_\mu, t_{\mu+1} the spline evaluation
    // depends on cpt_\mu-d..\cpt_\mu
    size_t sz = ks.size() - 1;
    cpts.assign(pts.cbegin() + sz - degree, pts.cbegin() + sz);
    std::copy(pts.cbegin(), pts.cbegin() + sz,std::back_inserter(cpts));


    // within the range, t_\mu, t_{\mu+1} the spline evaluation
    // depends on t_{\mu-d+1}.....t_\mu,t_{\mu+1},....t_{\mu+d}
    // (check the R_d matrix)

    // ./media/rmatrix.png

    // at \mu = 0, we need to ensure that t_{\mu-d+1} and
    // accessible, so d-1 extra knots are prefixed here at \mu =
    // ks.size()-2, we need to ensure that t_{\mu+d} is accessible,
    // so d-1 extra knots are appended here as well.
    std::vector<double> t;
    t.reserve( ks.size() + 2*degree );
    std::vector<double> startDiffs(degree+1), endDiffs(degree+1);

    // get the adjacent diffs between ks_0,ks_1,...
    std::adjacent_difference(ks.cbegin(), ks.cbegin() + degree + 1,
                             startDiffs.begin());
    // get the adjacent diffs between ks_{n-1},ks_{n-2}, ....
    std::adjacent_difference(ks.rbegin(), ks.rbegin() + degree + 1,
                             endDiffs.begin());


    startDiffs[0] = ks.back();  //change to ks_{n-1} in preparation
    //for the partial sums below
    endDiffs[0] = ks.front();   //change the start of backwardly
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
    // ks.front() so the
    // +1
    t.insert(t.end(),back_plus.begin() + 1, back_plus.end()); //already
    //copied
    //ks.back()
    return std::make_tuple(std::move(cpts),std::move(t),degree);
}

template <class CptsT, class KnotsT>
auto
make_periodic_bspline_wrap(
    const CptsT& pts,
    const KnotsT& ks,
    int degree_)
    ->periodic_bspline <RAWTYPE(pts[0]) >
{

    typedef RAWTYPE(pts[0]) point_t;
    return periodic_bspline < point_t >( std::move( make_bspline
                                                    (wrap_spline_params(pts, ks, degree_)) ));
}

template <class CptsT, class KnotsT>
auto
make_periodic_bspline (CptsT  pts,
                       KnotsT ks, int degree_)
    ->periodic_bspline <RAWTYPE(pts[0]) >
{
    typedef RAWTYPE(pts[0]) point_t;
    return periodic_bspline < point_t >( std::move( make_bspline (std::forward < CptsT > (pts),
                                                                  std::forward < KnotsT > (ks),
                                                                  degree_) ));
}

template <class Point>
periodic_bspline<Point>
make_periodic_bspline( bspline<Point>&& spl)
{
    assert(qry::is_periodic(spl));
    return periodic_bspline<Point>(std::forward<bspline<Point>>(spl));
}

}
#endif // ASTI_PERIODIC_BSPLINE_CONS

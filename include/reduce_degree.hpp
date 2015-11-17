#ifndef ASTI_REDUCE_DEGREE
#define ASTI_REDUCE_DEGREE
#include "split_into_bezier_patches.hpp"
#include "bspline_queries.hpp"
#include "reparametrize.hpp"
#include "change_basis.hpp"
#include "legendre_form.hpp"
#include <functional>
#include <type_traits>
#include <list>
#include "smat.hpp"
#include "trim_extend_join.hpp"
#include "bspline_x_cons.hpp"
#include "range.hpp"
#include "util.hpp"
#include <iterator>
#include "modify_bspline.hpp"
#include "reverse_curve.hpp"

namespace geom {
namespace ops {



// get a bezier spline of degee `deg' that is the closest L_2
// approximation of the given curve `crv'.  Here we use the nice
// properties of legendre polynomials for the L_2 approximation. Set
// the keep_* variables to constrain the resulting curve to pass
// through start or end
template <class SplineCurve>
SplineCurve
reduce_bezier_degree(const SplineCurve& crv, int deg,
                     bool keep_start = false,
                     bool keep_end = false)
{
    assert(qry::is_bezier(crv));
    auto const &  legf = to_legendre(qry::get_spline(crv));
    typedef RAWTYPE(legf) legf_t;
    auto cfs = legf.coeffs();
    cfs.erase(cfs.begin() + deg + 1, cfs.end());
    assert(cfs.size());
    return to_bezier(legf_t(cfs,legf.start_param(),
                            legf.end_param()));
}

// get a bspline curve of degree `deg' that is the closest (almost
// everywhere) L_2 approximation of the given curve `crv'.
template <class SplineCurve>
SplineCurve
reduce_degree(const SplineCurve& crv, int deg)
{
    assert(deg <= crv.degree());
    if(deg == crv.degree())
        return crv;
    auto const &patches =
        split_into_bezier_patches(qry::get_spline(crv));

    RAWTYPE(patches) sdpatches;
    for(auto const &p : patches) {
        sdpatches.push_back(reduce_bezier_degree(p,deg));
    }

    int i = 0;
    auto pf = sdpatches.begin();
    auto pe = sdpatches.end();

    auto nc = *pf;
    if(pf!=pe)
    {
        while(++pf!=pe)
        {
            double s = qry::start_param(nc);
            double e = qry::end_param(*pf);
            reparametrize(join_starts(reverse_curve(std::move(nc)),*pf,deg-1),s,e).swap(nc);
        }
    }

    typedef typename SplineCurve::cpts_t::iterator cref_t;
    typedef typename SplineCurve::knots_t::iterator kref_t;
    auto tweak_ends_reg = [&nc,&crv](  )
        {
            return SplineCurve(
                modify_bspline(nc,
                [&crv](cref_t pb, cref_t pe, kref_t, kref_t){
                    if(pb != pe) {
						*pb     = crv.eval(qry::start_param(crv));
						*(--pe) =  crv.eval(qry::end_param(crv));
                    }
                }) );
        };

    auto tweak_ends_per = [&nc,&deg,&crv](  )
        {
            modify_bspline(nc, [&crv](cref_t pb, cref_t pe, kref_t, kref_t) {
                    if(pb != pe) {
                        *pb     = crv.eval(qry::start_param(crv));
						*(--pe) =  crv.eval(qry::end_param(crv));
                    }
                });

            size_t sz = nc.knots().size();
            auto const &ks = nc.knots();
            std::vector<double> startDiffs(deg + 1), endDiffs(deg + 1);


            std::adjacent_difference(ks.cbegin() + deg, ks.cbegin() + 2*deg + 1,
                                     startDiffs.begin());

            std::adjacent_difference(ks.rbegin() + deg, ks.rbegin() + 2*deg + 1,
                                     endDiffs.begin());

            startDiffs[0] = ks.back();
            endDiffs[0] = ks.front();

            //  ./media/periodic.png

            std::vector<double> back_plus(startDiffs.size()),
            front_minus(endDiffs.size());

            std::partial_sum( startDiffs.begin(),
                              startDiffs.end(),
                              back_plus.begin() );

            std::partial_sum( endDiffs.begin(),
                              endDiffs.end(),
                              front_minus.begin() );

            rebase_at_start(nc, front_minus.rbegin()).swap(nc);
            rebase_at_end(nc,   back_plus.begin()).swap(nc);


            return SplineCurve( modify_bspline(nc,
                   [&deg](cref_t pb,
                           cref_t pe,
                           kref_t,
                           kref_t) {
                        auto p = pb;
                        for(auto q = pe - deg;
                            p != pb + deg;
                            ++p, ++q)
                        {
                            auto const &a = *p;
                            auto const &b = *q;
                            *p = *q = lerp(0.5,a,b);
                        }
                    }) );
        };
    typedef typename spline_traits<SplineCurve>::ptag ptag;
    return util::tag_switcher::eval(tweak_ends_per,tweak_ends_reg,
                              ptag());
}
}
}
#endif // ASTI_REDUCE_DEGREE

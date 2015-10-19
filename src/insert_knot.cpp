#include <iterator>
#include "rmat.hpp"
#include "smat.hpp"
#include "insert_knot.hpp"
#include "bspline_x_cons.hpp"
#include "util.hpp"
#include <numeric>
#include "type_utils.hpp"
#include "reparametrize.hpp"

namespace geom {

//{{{ --(@* "merge knots keeping duplicates within each list")
std::vector<double>
merge_knots(const std::vector<double> &t1,
            const std::vector<double> &t2)
{
    auto first1 = t1.cbegin();
    auto first2 = t2.cbegin();

    auto last1 = t1.cend();
    auto last2 = t2.cend();
    std::vector<double> t;
    t.reserve(t1.size() + t2.size() );
    while(first1 != last1 && first2 != last2)
    {
        if(tol::eq(*first1,*first2))
        {
            t.push_back(*first1);
            first1++;
            first2++;
        }
        else if(*first1 < *first2)
        {
            t.push_back(*first1++);
        }else
            t.push_back(*first2++);
    }
    auto dest = std::copy(first1, last1, std::back_inserter(t));    // copy any tail
    std::copy(first2, last2, dest);
    return t;
}
//}}}

//{{{ --(@* "implements the oslo algorithm to insert knots")
// see(@url :file-name "http://www.uio.no/studier/emner/matnat/ifi/INF-MAT5340/v09/undervisningsmateriale/book.pdf#page=100" :display "page100")

template <class SplineType>
SplineType
ops::insert_knots(const SplineType& crv,
                  const std::vector<double>& new_knots)
{
    typedef typename SplineType::knots_t knots_t;

    knots_t taus;
    auto const & t = crv.knots();

    taus.reserve( new_knots.size() + t.size() );
    auto const & cpts = crv.control_points();
    typedef decltype(cpts) cpts_t;
    taus = merge_knots(t,new_knots);
    if(taus.size() == t.size() )
        return SplineType(crv);

    typedef RAWTYPE(cpts[0]) point_t;
    rmat<point_t> m(cpts, t, crv.degree());
    cpts_t new_cpts(m.insert_knots(taus));
    return make_bspline(
        std::move(new_cpts),
        std::move(taus),
        crv.degree(),
        typename spline_traits<SplineType>::ptag(),
        typename spline_traits<SplineType>::rtag())
        ;
}

// ./media/boehm.png
template <class SplineType>
SplineType
insert_knot_impl( const SplineType &crv,
                  double u,
                  geom::regular_tag )
{
    typedef typename SplineType::point_t point_t;
    typedef typename SplineType::knots_t knots_t;
    auto & cpts = crv.control_points();
    auto p = crv.degree();
    auto & t = crv.knots();
    size_t nu = rmat_base_vd(t,p).locate_nu(u);
    typename SplineType::cpts_t newcpts;
    newcpts.reserve(cpts.size() + 1);
    std::copy(cpts.cbegin(), cpts.cbegin() + nu - p + 1,
              std::back_inserter(newcpts));
    for(size_t j = nu - p + 1;j < nu + 1; ++j)
    {
        auto lambda = sdiv(t[j + p] - u, t[j + p] - t[j]);
        newcpts.push_back( lerp(lambda, cpts[j], cpts[j - 1]) );
    }
    std::copy(cpts.cbegin() + nu, cpts.cend(),
              std::back_inserter(newcpts));
    typename SplineType::knots_t newts;
    newts.reserve(t.size() + 1);
    newts.assign(t.cbegin(),t.cbegin() + nu + 1);
    newts.push_back(u);
    newts.insert(newts.end(), t.cbegin() + nu + 1, t.cend());
    return make_bspline (
        std::move(newcpts),std::move(newts),p,
        typename spline_traits<SplineType>::ptag(),
        typename spline_traits<SplineType>::rtag())
        ;
}


template <class SplineType>
SplineType
insert_knot_impl( const SplineType &crv,
                  double u,
                  geom::periodic_tag)
{
    u               = periodic_param(crv.param_range(), u);
    auto         p  = crv.degree();
    auto &       t  = crv.knots();
    size_t       nu = rmat_base_vd(t,p).locate_nu(u);

    auto crv0 = insert_knot_impl(crv, u, geom::regular_tag());
    if( nu < size_t(2*p) )
    {
        // insertion must be symmetric
        std::vector<double> startDiffs(p+1);
        auto &t0  = crv0.knots();
        std::adjacent_difference(t0.cbegin() + p, t0.cbegin()
                                 + 2*p  + 1,
                                 startDiffs.begin());
        startDiffs[0] = t0.cend()[-p-1];
        std::vector<double> ks(p+1);
        std::partial_sum(startDiffs.cbegin(),
                         startDiffs.cend(), ks.begin());
        return rebase_at_end(crv0, ks.cbegin());
    }else if( nu > t.size() - 2 * p - 1)
    {
        std::vector<double> endDiffs(p+1);
        auto &t0 = crv0.knots();
        std::adjacent_difference(t0.cend() - 2*p - 1,
                                 t.cend() - p,
                                 endDiffs.begin());
        endDiffs[0] = t0[p];
        std::vector<double> ks(p+1);
        std::partial_sum(endDiffs.cbegin(),
                         endDiffs.cend(), ks.begin());
        return rebase_at_start(crv0, ks.cbegin());
    }
    return crv0;
}

template <class SplineType>
SplineType
ops::insert_knot(const SplineType& crv,
                 double u)
{
    return insert_knot_impl(crv,u,
                            typename spline_traits<SplineType>::ptag());
}

//}}}


//{{{
template <class SplineType>
bool
ops::match_knots(SplineType& spl1, SplineType& spl2)
{
    assert(spl1.degree() == spl2.degree());
    auto const & rspl1 = reparametrize(spl1); // get a 0, 1 parametrization
    auto const & rspl2 = reparametrize(spl2);

    insert_knots(rspl1, rspl2.knots()).swap(spl1);
    insert_knots(rspl2, rspl1.knots()).swap(spl2);
    return true;
}

//}}}
}
/*
  Local Variables:
  eval:(load-file "./scripts/temp.el")
  eval:(setq methods (list "insert_knots"
  "insert_knot" "match_knots"
  ))
  eval:(setq spltypes (list "bspline<double>"
  "bspline<point2d_t>"
  "bspline<point3d_t>"
  "bspline<point4d_t>"
  "periodic_bspline<double>"
  "periodic_bspline<point2d_t>"
  "periodic_bspline<point3d_t>"
  "periodic_bspline<point4d_t>"
  "rational_bspline < point2d_t,regular_tag>"
  "rational_bspline < point3d_t,regular_tag>"
  "rational_bspline < double, regular_tag>"
  "rational_bspline < point2d_t,periodic_tag>"
  "rational_bspline < point3d_t,periodic_tag>"
  "rational_bspline < double,periodic_tag>"
  ))
  eval:(instantiate-templates "insert_knot" "ops" (list ) (product
  methods spltypes ))
  End:
  // dump all explicitly instantiated templates below
  */
#include "bspline.hpp"
#include "periodic_bspline.hpp"
#include "rational_bspline.hpp"
#include "point.hpp"
namespace geom {
#include "insert_knot_inst.inl"
}

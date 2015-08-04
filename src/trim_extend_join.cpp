//-*- mode:c++ -*-
#include "point.hpp"
#include "trim_extend_join.hpp"
#include "constant_iterator.hpp"
#include "raise_degree.hpp"
#include "tol.hpp"
#include "reparametrize.hpp"
#include "geom_exception.hpp"
#include "reverse_curve.hpp"
#include "insert_knot.hpp"
#include "smat.hpp"

namespace geom {

template <class SplineType>
SplineType
ops::join_starts(const SplineType& spl1,
                 const SplineType& spl2,
                 int join_cont)
{
    SplineType spl1_clamped (clamp_end(spl1));
    SplineType spl2_clamped (clamp_start(spl2));

    // let the two splines have the same start parameters
    reparametrize_start(spl1_clamped, 0).swap(spl1_clamped);
    reparametrize_start(spl2_clamped, 0).swap(spl2_clamped);

    match_degrees(spl1_clamped, spl2_clamped);

    int p = spl1_clamped.degree();
    if(join_cont >= p )
        throw geom_exception(continuity_condition_too_tight);

    typename SplineType::knots_t ks(p + 1);

    auto a = spl1_clamped.knots().begin()[p+1];
    auto b = spl2_clamped.knots()[p+1];

    if(a > b ) std::swap(a,b);
    std::fill_n(ks.begin(), join_cont + 1, -a);
    std::fill_n(ks.begin() + join_cont + 1, p - join_cont,0);

    rmat_base_vd r(spl2_clamped.knots(),p);
    ptrdiff_t nk = join_cont + 1 - r.mult(a);
    for( ;nk > 0; --nk ) {
        insert_knot(spl2_clamped,a).swap(spl2_clamped);
    }
    rmat_base_vd(spl1_clamped.knots(),p).swap(r);
    nk = join_cont + 1 - r.mult(a);
    for( ;nk > 0; --nk ) {
        insert_knot(spl1_clamped,a).swap(spl1_clamped);
    }

    rebase_at_start(spl1_clamped, ks.begin()).swap(spl1_clamped);
    rebase_at_start(spl2_clamped, ks.begin()).swap(spl2_clamped);

    auto & cpts1 = spl1_clamped.control_points();
    auto & cpts2 = spl2_clamped.control_points();
    typename SplineType::cpts_t cpts;
    typedef typename SplineType::cpts_t::const_reference cref;
    typedef typename SplineType::point_t point_t;
    cpts.reserve(cpts1.size() + cpts2.size());
    std::copy(cpts1.crbegin(), cpts1.crend(), std::back_inserter(cpts));
    // merge the p + 1 cpts at the end of reversed c0 and start of c1
    std::transform(cpts.cend() - (join_cont + 1), cpts.cend(),
                   cpts2.cbegin(),/* into */ cpts.end() - (join_cont + 1),
                   [](cref p1,cref p2){
                       return lerp(0.5, p1,p2);} );

    std::copy(cpts2.cbegin() + join_cont + 1,cpts2.cend(),
              std::back_inserter(cpts));
    // merge the two knot sequences together to get the knot sequence
    // for the join
    std::vector<double> newknots;
    newknots.reserve(spl1_clamped.knots().size() +
                     spl2_clamped.knots().size() - p - 1);

    reverse_curve(spl1_clamped).swap(spl1_clamped);
    newknots.assign(spl1_clamped.knots().cbegin(),
                    spl1_clamped.knots().cend() - (join_cont + 1));
    newknots.insert(newknots.end(),
                    spl2_clamped.knots().cbegin() + (p + 1),
                    spl2_clamped.knots().cend());

    typedef spline_traits<SplineType> str;
    return make_bspline (std::move(cpts),
                         std::move(newknots), p,
                         typename str::ptag(),
                         typename str::rtag()
        );
}


template <class SplineType>
SplineType
ops::extract_regular_curve(const SplineType &spl)
{
    auto &ts    = spl.knots();
    int   deg   = spl.degree();

    assert(!tol::param_eq(ts.front(), ts.back()));
    return clamp_end(clamp_start(spl));
}


template <class SplineType>
SplineType
ops::trim_curve(const SplineType &spl, double a, double b)
{
    auto &ts    = spl.knots();
    int   deg   = spl.degree();

    assert(!tol::param_eq(ts.front(), ts.back()));
    return clamp_at_right(b, clamp_at_left(a, spl));
}

template <class SplineType>
SplineType
ops::extend_curve_start(const SplineType & spl, double delta)
{
    auto const & s = clamp_start(spl);
    auto pr =  s.param_range();
    double u = pr.first;
    u -= delta * (pr.second - pr.first);
    return rebase_at_start(s, util::make_constant_iterator(u) );
}

template <class SplineType>
SplineType
ops::extend_curve_end(const SplineType & spl, double delta)
{
    auto const & s  = clamp_start(spl);
    auto         pr = s.param_range();
    double       v  = pr.second;
    v += (delta * (pr.second -  pr.first));
    return rebase_at_end(s, util::make_constant_iterator(v));
}

template <class SplineType>
SplineType
ops::extend_curve_end_to_pt(const SplineType & spl,
                            typename SplineType::point_t const & target)
{
    auto const & s         = reparametrize (clamp_start(spl), 0, 1);
    auto const & t         = s.knots();
    size_t       n         = s.control_points().size();
    int          d         = s.degree();
    double       chord_len = 0;
    auto         pt         = s.eval(t[d]);

    for(size_t r = 0;r < n - d; ++r) {
        auto newc = s.eval(t[d + r + 1]);
        chord_len += len(newc - pt);
        pt = newc;
    }

    double ld  = len(target - pt);
    chord_len += ld;
    auto delta      = ld / chord_len;

    std::vector<double> ks(d + 1, 1 + delta);
    ks[0] = 1;

    auto exs(rebase_at_end(s, ks.begin()));
    auto newks(exs.knots());
    auto newcpts(exs.control_points());
    newks.push_back(1 + delta);
    typedef typename SplineType::cpts_t::value_type cpt_val_t;
    newcpts.push_back(cpt_val_t(target));

    typedef spline_traits<SplineType> str;
    return make_bspline
        (std::move(newcpts), std::move(newks), d,
         typename str::ptag(), typename str::rtag());

}



}
//{{{  instantiation scripts

/*
  Local Variables:
  eval:(load-file "./scripts/temp.el")
  eval:(setq methods (list "extract_regular_curve"
  "trim_curve"
  "extend_curve_start"
  "extend_curve_end"
  "join_starts"
  "extend_curve_end_to_pt"
  ))
  eval:(setq spltypes (list "bspline<double>"
  "bspline<point2d_t>"
  "bspline<point3d_t>"
  "bspline<point4d_t>"
  "rational_bspline < point2d_t,regular_tag>"
  "rational_bspline < point3d_t,regular_tag>"
  "rational_bspline < double, regular_tag>"
  ))
  eval:(instantiate-templates "trim_extend_join" "ops" (list )
  (product methods spltypes ))
  End:
  // dump all explicitly instantiated templates below
  */
//}}}

//{{{  instantiation
#include "bspline.hpp"
#include "rational_bspline.hpp"
#include "point.hpp"
namespace geom {
#include "trim_extend_join_inst.inl"
}
//}}}

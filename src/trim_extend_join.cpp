#include "trim_extend_join.hpp"
#include "constant_iterator.hpp"
#include "smat.hpp"
#include "raise_degree.hpp"
#include "tol.hpp"
#include "reparametrize.hpp"
#include "geom_exception.hpp"
#include "reverse_curve.hpp"
#include "insert_knot.hpp"
namespace geom {

template <class SplineType>
SplineType
bspline_ops::join_starts(const SplineType& spl1,
                         const SplineType& spl2,
                         int join_cont)
{
    // assert(!spl1.is_periodic() && !spl2.is_periodic());
    SplineType spl1_clamped (clamp_end(spl1));
    SplineType spl2_clamped (clamp_start(spl2));

    // let the two splines have the same start parameters
    spl1_clamped.swap(reparametrize_start(spl1_clamped, 0));
    spl2_clamped.swap(reparametrize_start(spl2_clamped, 0));

    int p1 = spl1_clamped.degree();
    int p2 = spl2_clamped.degree();
    for(;p2 < p1;++p2)
        raise_degree(spl2_clamped);

    for(;p1 < p2;++p1)
        raise_degree(spl1_clamped);

    int p = spl1_clamped.degree();
    if(join_cont >= p )
        throw geom_exception(continuity_condition_too_tight);
    SplineType::knots_t ks(p + 1);

    auto a = spl1_clamped.knots().begin()[p+1];
    auto b = spl2_clamped.knots()[p+1];
    if(a > b ) std::swap(a,b);
    std::fill_n(ks.begin(), join_cont + 1, -a);
    std::fill_n(ks.begin() + join_cont + 1, p - join_cont,0);
    
    // match the widths of the first knot interval to 'a'
    //spl2_clamped.swap(reparametrize(spl2_clamped, 0, a / b));

    rmat_base_vd r(spl2_clamped.knots(),p);
    ptrdiff_t nk = join_cont + 1 - r.mult(a);
    for( ;nk > 0; --nk ) {
        spl2_clamped.swap(insert_knot(spl2_clamped,a));
    }
    r.swap(rmat_base_vd(spl1_clamped.knots(),p));
    nk = join_cont + 1 - r.mult(a);
    for( ;nk > 0; --nk ) {
        spl1_clamped.swap(insert_knot(spl1_clamped,a));
    }

    spl1_clamped.swap(rebase_at_start(spl1_clamped, ks.begin()));
    spl2_clamped.swap(rebase_at_start(spl2_clamped, ks.begin()));

    auto & cpts1 = spl1_clamped.control_points();
    auto & cpts2 = spl2_clamped.control_points();
    typename SplineType::cpts_t cpts;
    typedef typename SplineType::cpts_t::const_reference cref;
    typedef SplineType::point_t point_t;
    cpts.reserve(cpts1.size() + cpts2.size());
    std::copy(cpts1.crbegin(), cpts1.crend(), std::back_inserter(cpts));
    // merge the p + 1 cpts at the end of reversed c0 and start of c1
    std::transform(cpts.cend() - (join_cont + 1) ,
                   cpts.cend(),
                   cpts2.cbegin(),
                   cpts.end() - (join_cont + 1),
                   [](cref p1,cref p2){
                       return lerp(0.5, p1,p2);} );    

    std::copy(cpts2.cbegin() + join_cont + 1,cpts2.cend(),
              std::back_inserter(cpts));
    // merge the two knot sequences together to get the knot sequence
    // for the join
    std::vector<double> newknots;
    newknots.reserve(spl1_clamped.knots().size() +
                     spl2_clamped.knots().size() - p - 1);

    spl1_clamped.swap(reverse_curve(spl1_clamped) );

    newknots.assign(spl1_clamped.knots().cbegin(),
        spl1_clamped.knots().cend() - (join_cont + 1));
    newknots.insert(newknots.end(),
        spl2_clamped.knots().cbegin() + (p + 1) ,
                    spl2_clamped.knots().cend());
   
    return SplineType(std::move(cpts), std::move(newknots), p);
}


template <class SplineType>
SplineType
bspline_ops::extract_regular_curve(const SplineType &spl)
{
    auto &ts    = spl.knots();
    int   deg   = spl.degree();

    assert(!tol::param_eq(ts.front(), ts.back()));
    return clamp_end(clamp_start(spl));
}

template <class SplineType>
SplineType
bspline_ops::trim_curve(const SplineType &spl, double a, double b)
{
    auto &ts    = spl.knots();
    int   deg   = spl.degree();

    assert(!tol::param_eq(ts.front(), ts.back()));
    return clamp_at_right(b, clamp_at_left(a, spl));
}

template <class SplineType>
SplineType
bspline_ops::extend_curve_start(const SplineType & spl,
                                double delta)
{
    auto const & s = clamp_start(spl);
    double u = s.param_range().first;
    u -= delta;
    return rebase_at_start(s, util::make_constant_iterator(u) );
}

template <class SplineType>
SplineType
bspline_ops::extend_curve_end(const SplineType & spl,
                              double delta)
{
    auto const & s = clamp_start(spl);
    double v = s.param_range().second;
    v += delta;
    return rebase_at_end(s, util::make_constant_iterator(v));
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
  ))
  eval:(setq spltypes (list "bspline<double>"
  "bspline<point2d_t>"
  "bspline<point3d_t>"
  "bspline<point4d_t>"
  ))
  eval:(instantiate-templates "trim_extend_join" "bspline_ops" (list ) methods spltypes )
  End:
// dump all explicitly instantiated templates below
*/
//}}}

//{{{  instantiation
#include "bspline.hpp"
#include "periodic_bspline.hpp"
#include "point.hpp"
namespace geom {
#include "trim_extend_join_inst.cpp"
}
//}}}

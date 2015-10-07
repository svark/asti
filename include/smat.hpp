//-*- mode:c++ -*-
#ifndef ASTI_SMAT_HPP
#define ASTI_SMAT_HPP
#include <vector>
#include <algorithm>
#include "rmat.hpp"
#include "bspline_x_cons.hpp"
#ifndef NDEBUG
#include "reverse_curve.hpp"
#endif
//#include "util.hpp" // for constant_iterator
#include "constant_iterator.hpp"
namespace geom {
// bspline to bezier (and vice versa)  conversion  matrix.  as described in
// (@file :file-name "b2b.pdf" :to "b2b.pdf" :display "bspline2bezier")
//
// http://amsacta.unibo.it/853/1/preprint.pdf
//  bspline to bezier conversion and vice versa
struct smat : rmat_base_vd {
    double a; // interval start from which to extract curve, end will be
              // the next knot
    double b;
    smat(double a_, double b_, const knots_t& t_,
         int deg_) :rmat_base_vd(t_, deg_), a(a_),b(b_)
    {
    }

    template <class ArrayIter>
    void
    reval(ArrayIter cachea // expected to be of size degree + 1
        ) const
    {
        size_t nu = (locate_nu(a));
        int size = deg;
        auto t_ = t + nu;
        assert(nu >= size_t(size));

        typedef typename std::iterator_traits<ArrayIter>::value_type U;
        typedef decltype(mk_stdvec(U())) cpts_t;
        cpts_t cachec(size + 1);
        // todo:parallelize
#pragma loop(hint_parallel(8))
        for(int j = 0; j <= size; ++j) {
            cpts_t cacheb(cachea, cachea + size + 1);
            for(int sz = size; sz > j; --sz) {
                // ../src/media/reval.png
                double lambda = sdiv(t_[j + 1 - sz] - a, b - a);
                if(tol::param_eq(lambda,  0)) continue;
                for(int i = 0; i < sz; ++i){
                    cacheb[i] = lerp(lambda, cacheb[i], cacheb[i + 1]);
                }
            }

            for(int sz = j; sz >= 1; --sz) {
                // ../src/media/reval2.png
                double mu =  sdiv(t_[sz] - a, b - a);
                if(tol::param_eq(mu,  0)) continue;
                for(int i = 0; i < sz; ++i){
                    cacheb[i] = lerp(mu, cacheb[i], cacheb[i + 1]);
                }
            }
            cachec[j]=(cacheb[0]);
        }
        std::copy(cachec.begin(),cachec.end(),cachea);
    }

    template <class ArrayIter>
    void
    seval(ArrayIter cachea // expected to be of size: degree + 1
        ) const
    {
        size_t nu = (locate_nu(a));
        const int size = deg;
        auto t_ = t + nu;
        assert(nu >= size_t(size));
        // ./media/seval0.png
        typedef typename std::iterator_traits<ArrayIter>::value_type U;
        typedef decltype(mk_stdvec(U())) cpts_t;
        cpts_t cachec(size+1);
#pragma loop(hint_parallel(8))
        for(int i = 0; i <= size; ++i) {
            cpts_t cacheb(cachea, cachea + size + 1);
            for(int sz = size; sz > i; --sz) {
                for(int j = 1; j <= sz; ++j){
                    double lambda = sdiv(a - t_[j - sz],
                                         t_[j] - t_[j - sz]);
                    cacheb[j - 1] = lerp(lambda, (cacheb[j - 1]), (cacheb[j]));
                }
            }
            // ./media/seval1.png
            for(int sz = i; sz >= 1; --sz) {
                for(int j = 1; j <= sz; ++j){
                    double mu =  sdiv(b - t_[j - sz],
                                      t_[j] - t_[j - sz]);
                    cacheb[j - 1] = lerp(mu, (cacheb[j - 1]), (cacheb[j]));
                }
            }
            cachec[i] = cacheb[0];
        }
        std::copy(cachec.begin(), cachec.end(), cachea);
    }
};


template <class SplineCurve, class KnotIter>
SplineCurve rebase_at_left(const SplineCurve & crv,
                           double a,
                           KnotIter us)
{

    auto const & t = crv.knots();
    int deg  = crv.degree();
    auto const & cpts = crv.control_points();
    rmat_base_vd r(t, deg);
    size_t nu = r.locate_nu(a);
    // all the knots in us are expected to be <= a
#ifndef NDEBUG
    std::for_each(us, us + (deg + 1), [&a](double u){ assert(u <= a);});
#endif
    typedef typename SplineCurve::cpts_t cpts_t;
    typedef typename SplineCurve::knots_t knots_t;
    double b = t[nu+1];
    cpts_t newcpts(cpts.cbegin() + (nu - deg), cpts.cend());
    // get control points wrt bernstein basis
    smat(a, b, t, deg).seval(newcpts.begin());
    // switch back to bspline basis
    knots_t ks(t.cbegin() + (nu - deg), t.cend());
    std::copy_n(us, deg + 1, ks.begin()); // knots updated to reclamp at
    // 'a'
    smat(a, b, ks, deg).reval(newcpts.begin());

    typedef spline_traits<SplineCurve> str;
    return make_bspline
        (std::move(newcpts), std::move(ks), deg,
         typename str::ptag(),
         typename str::rtag())
        ;
}

template <class SplineCurve, class FnType>
SplineCurve transform_at_left(const SplineCurve & crv,
                              FnType f /*std::function<void(PointIter&)>*/
    )
{
    double a = crv.param_range().first;
    int deg  = crv.degree();
    auto const & cpts = crv.control_points();
    double b = crv.knots()[deg+1];
    transform_at(crv,a,b,f);
}

// modify curve by changing the bezier control points between a,b of
// given curve. New curve has the same knots as the old.
template <class SplineCurve, class FnType>
SplineCurve transform_at(const SplineCurve & crv,
                         double a, double b,
                         FnType f /*std::function<void(PointIter&)>*/
    )
{
    auto const & t = crv.knots();
    int deg  = crv.degree();
    auto const & cpts = crv.control_points();
    rmat_base_vd r(t, deg);
    size_t nu = r.locate_nu(a);

    typedef typename SplineCurve::cpts_t cpts_t;
    typedef typename SplineCurve::knots_t knots_t;

    cpts_t newcpts(cpts);
    // get control points wrt bernstein basis
    smat(a, b, t, deg).seval(newcpts.begin());
    // switch back to bspline basis after transforming the points
    f(newcpts.begin(),a,b);
    smat(a, b, t, deg).reval(newcpts.begin());

    typedef spline_traits<SplineCurve> str;
    return make_bspline
        (std::move(newcpts), t, deg,
         typename str::ptag(),
         typename str::rtag())
        ;
}

template <class SplineCurve, class KnotIter>
SplineCurve rebase_at_right(const SplineCurve & crv,
                            double b, KnotIter us)
{

    auto & t = crv.knots();
    int deg  = crv.degree();
    auto const & cpts = crv.control_points();
    rmat_base_vd r(t, deg);

    size_t nu = r.locate_nu(b - tol::param_tol/2);
    //  assert(r.locate_nu(b+tol::param_tol/2)==nu+1);
    double a = t[nu];
    // all the knots in us are expected to be >= b
#ifndef NDEBUG
    std::for_each(us, us + (deg + 1), [&b](double u){ assert(u >= b);});
#endif

    typedef typename SplineCurve::cpts_t cpts_t;
    typedef typename SplineCurve::knots_t knots_t;
    cpts_t newcpts(cpts.cbegin(), cpts.cbegin() + (nu + 1));
    // get deg + 1, control points wrt bernstein basis
    smat(a, b, t, deg).seval(newcpts.begin() + (nu - deg));
    // switch back to bspline basis
    knots_t ks(t.cbegin(), t.cbegin() + (nu + deg + 2));
    std::copy_backward(us, us + (deg + 1), ks.end());
    smat(a, b, ks, deg).reval(newcpts.begin() + (nu - deg));

    typedef spline_traits<SplineCurve> str;
    return make_bspline (std::move(newcpts),
                         std::move(ks), deg,
                         typename str::ptag(),
                         typename str::rtag())
        ;
}

template <class SplineCurve, class KnotIter>
SplineCurve rebase_at_start(const SplineCurve & crv, KnotIter us)
{

    auto & t = crv.knots();
    int deg  = crv.degree();
    rmat_base_vd r(t, deg);
    double a = crv.param_range().first;
    return rebase_at_left(crv, a, us);
}

template <class SplineCurve, class KnotIter>
SplineCurve rebase_at_end(const SplineCurve & crv, KnotIter us)
{

    auto & t = crv.knots();
    int deg  = crv.degree();
    rmat_base_vd r(t, deg);
    double b = crv.param_range().second;
    return rebase_at_right(crv, b, us);
}

template <class SplineCurve>
SplineCurve clamp_at_right(double b, const SplineCurve & crv)
{
    return rebase_at_right(crv, b, util::make_constant_iterator(b));
}

template <class SplineCurve>
SplineCurve clamp_at_left(double a, const SplineCurve & crv)
{
    return rebase_at_left(crv, a, util::make_constant_iterator(a));
}

template <class SplineCurve>
SplineCurve clamp_start(const SplineCurve & crv)
{
    auto & t = crv.knots();
    int deg  = crv.degree();
    rmat_base_vd r(t, deg);

    if(r.start_mult() == deg + 1)
        return SplineCurve(crv);

    double a = t[deg];
    return clamp_at_left(a,crv);
}

template <class SplineCurve>
SplineCurve clamp_end(const SplineCurve & crv)
{
    auto & t = crv.knots();
    int deg  = crv.degree();
    rmat_base_vd r(t, deg);

    if(r.end_mult() == deg + 1)
        return SplineCurve(crv);

    double b = t[t.size() - deg - 1];
    return clamp_at_right(b, crv);
}

#ifndef NDEBUG
// for testing
// yaim -> yet another impl

template <class SplineCurve>
SplineCurve clamp_at_right_yaim(double b, const SplineCurve & crv)
{
    typedef typename SplineCurve::cpts_t cpts_t;
    SplineCurve rcrv(ops::reverse_curve(crv));
    return ops::reverse_curve(clamp_at_left(-b, rcrv));
}

template <class SplineCurve>
SplineCurve clamp_end_yaim(const SplineCurve & crv)
{
    double b = crv.param_range().second;
    return clamp_at_right_yaim(b, crv);
}

template <class SplineCurve>
SplineCurve clamp_at_left_yaim(double a, const SplineCurve & crv)
{
    typedef typename SplineCurve::cpts_t cpts_t;
    SplineCurve rcrv(ops::reverse_curve(crv));
    return ops::reverse_curve(clamp_at_right(-a, rcrv));
}

template <class SplineCurve>
SplineCurve clamp_start_yaim(const SplineCurve & crv)
{
    double a = crv.param_range().first;
    return clamp_at_left_yaim(a, crv);

}
#endif

}
#endif //SMAT_HPP

//-*- mode:c++ -*-
#ifndef SMAT_HPP
#define SMAT_HPP
#include <vector>
#include <list>
#include "rmat.hpp"

namespace geom {
// bspline to bezier (and vice versa)  conversion  matrix.  as described in
// (@file :file-name "b2b.pdf" :to "b2b.pdf" :display "bspline2bezier")
//

struct  bspline_2_bezier_mat_t{};
struct  bezier_2_bspline_mat_t{};

//  bspline to bezier conversion and vice versa
struct smat : rmat_base_vd{
    double a; // interval start from which to extract curve, end will
              // be the next knot

    smat(double a_, const knots_t& t_,
         int deg_) :rmat_base_vd(t_, deg_), a(a_)
    {
    }

    template <class T>
    void
    reval(T cachea // expected to be of size degree + 1
        ) const
    {
        double u = a;
        size_t nu = (locate_nu(a));
        double b = t[nu + 1];
        int size = deg;
        auto t_ = t + nu;
        assert(nu >= size);

        typedef typename std::decay < decltype(cachea[0]) >::type U;
        std::vector<U> cachec(size + 1);
        // todo:parallelize
        for(int j = 0; j <= size; ++j) {
            std::vector<U> cacheb(cachea, cachea + size + 1);
            for(int sz = size; sz > j; --sz) {
       // ./reval.png
                double lambda = sdiv(t_[j + 1 - sz] - a, b - a);
                if(tol::param_eq(lambda,  0)) continue;
                for(int i = 0; i < sz; ++i){
                    cacheb[i] = lerp(lambda, cacheb[i], cacheb[i + 1]);
                }
            }

            for(int sz = j; sz >= 1; --sz) {
       // ./reval2.png
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

    template <class T>
    void
    seval(T cachea // expected to be of size: degree + 1
         ) const
    {
        double u = a;
        size_t nu = (locate_nu(a));
        double b = t[nu + 1];
        int size = deg;
        auto t_ = t + nu;
        assert(nu >= size);
        // ./seval0.png
        typedef typename std::decay < decltype(cachea[0]) >::type U;
        std::vector<U> cachec(size+1);
        for(int i = 0; i <= size; ++i) {
            std::vector<U> cacheb(cachea, cachea + size + 1);
            for(int sz = size; sz > i; --sz) {
                for(int j = 1; j <= sz; ++j){
                    double lambda = sdiv(a - t_[j - sz],
                                         t_[j] - t_[j - sz]);
                    cacheb[j - 1] = lerp(lambda, cacheb[j - 1], cacheb[j]);
                }
            }
            // ./seval1.png
            for(int sz = i; sz >= 1; --sz) {
                for(int j = 1; j <= sz; ++j){
                    double mu =  sdiv(b - t_[j - sz],
                                  t_[j] - t_[j - sz]);
                    cacheb[j - 1] = lerp(mu, cacheb[j - 1], cacheb[j]);
                }
            }
            cachec[i] = cacheb[0];
        }
        std::copy(cachec.begin(), cachec.end(), cachea);
    }
};

template <class SplineCurve>
SplineCurve clamp_at_left(double a, const SplineCurve & crv)
{
    auto & t = crv.knots();
    int deg  = crv.degree();
    auto const & cpts = crv.control_points();
    rmat_base_vd r(t, deg);
    size_t nu = r.locate_nu(a);

    typedef SplineCurve::cpts_t cpts_t;
    typedef SplineCurve::knots_t knots_t;

    cpts_t newcpts(cpts.cbegin() + nu - deg, cpts.cend());
    // get control points wrt bernstein basis
    smat(a, t, deg).seval(newcpts.begin());
    // switch back to bspline basis
    knots_t ks(t.cbegin() + nu - deg, t.cend());
    std::fill_n(ks.begin(), deg + 1, a); // knots updated to clamp at
                                         // 'a'
    smat(a, ks,
         deg).reval(newcpts.begin());

    return SplineCurve(std::move(newcpts), std::move(ks), deg)
        .translate(crv.base_point());
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
SplineCurve clamp_at_right(double b, const SplineCurve & crv)
{
    typedef SplineCurve::cpts_t cpts_t;
    SplineCurve rcrv(bspline_ops::reverse_curve(crv));
    return bspline_ops::reverse_curve(clamp_at_left(-b, rcrv));
}

template <class SplineCurve>
SplineCurve clamp_end(const SplineCurve & crv)
{
    auto & t = crv.knots();
    int deg  = crv.degree();
    rmat_base_vd r(t, deg);

    if(r.end_mult() == deg + 1)
        return SplineCurve(crv);
    typedef SplineCurve::cpts_t cpts_t;
    SplineCurve rcrv(bspline_ops::reverse_curve(crv));
    return bspline_ops::reverse_curve(clamp_start(rcrv));
}

}
#endif //SMAT_HPP

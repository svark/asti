//-*- mode:c++ -*-
#ifndef ASTI_CHANGE_BASIS
#define ASTI_CHANGE_BASIS
#include "point.hpp"
#include "bspline.hpp"
#include "bspline_queries.hpp"
#include "type_utils.hpp"
#include "bezier_form.hpp"
#include "monomial_form.hpp"
#include "legendre_form.hpp"

namespace geom {
namespace ops  {

template<class Point>
monomial_form < Point >
to_monomial(const bspline < Point >& bezf)
{
    assert(qry::is_bezier(bezf));
    auto & b = bezf.control_points();
    size_t sz =  b.size();
    typedef  typename monomial_form < Point >::cpts_t cpts_t;
    cpts_t monf(sz, Point(0.0));
    long cni = 1;
    short signi = 1;
    int n = sz - 1;
    for(size_t i =  0;i < sz; ++i, signi = -signi) {
        if(i != 0){
            cni *= (n - i + 1) ;
            cni /= i;
        }
        long cil =  1;
        short signl = 1;
        for(size_t l =  0;l <= i; ++l, signl = -signl)
        {
            if(l != 0) {
                cil *= (i - l + 1);
                cil /= l;
            }
            auto const & ev(b[l]);
            double co =  signi * signl * cni * cil ;
            monf[i] += (scaled_copy(ev, co));
        }
    }
    double s,e;
    std::tie(s,e) = bezf.param_range();
    return monomial_form < Point > (std::move(monf), s, e);
}

inline
std::vector<size_t> ncks(size_t n)
{
    std::vector<size_t> ncks_;
    ncks_.reserve(n + 1);
    size_t val = 1;
    ncks_.push_back(1);
    for(size_t k = 1; k <= n; ++k)
    {
        val *= (n - k + 1);
        val /= k;
        ncks_.push_back(val);
    }
    return ncks_;
}

inline
std::vector<size_t> nkcks(size_t n, size_t j)
{
    std::vector<size_t> nkcks_;
    nkcks_.reserve(j + 1);
    size_t val = 1;
    nkcks_.push_back(1);
    for(size_t k = 1; k <= j; ++k)
    {
        val *= (n + k);
        val /= k;
        nkcks_.push_back(val);
    }
    return nkcks_;
}

template<class Point>
bezier_form<Point>
to_bezier(const monomial_form < Point > & mf )
{
    int n  =  mf.size() - 1;
    typename bezier_form<Point>::cpts_t  cpts(n + 1, Point(0.0));
    for(size_t l = 0;l < mf.size(); ++l) {
        long clk = 1, cnk = 1;
        for(size_t k = 0; k <= l; ++k)
        {
            if(k != 0) {
                cnk *= (n - k + 1) ;
                cnk /= k;
                clk *= (l - k + 1);
                clk /= k;
            }
            cpts[l] += double(clk) * mf[k]/ cnk;
        }
    }
    return make_bezier_form(std::move(cpts),
                            mf.start_param(), mf.end_param());
}
    // ../src/media/to_legendre.png
    template <class Point>
    legendre_form<Point>
    to_legendre(const bspline<Point>& bezf)
    {
        assert(qry::is_bezier(bezf));
        auto & b   = bezf.control_points();
        size_t sz  = b.size();
        auto &cpts = bezf.control_points();
        auto &ks   = bezf.knots();
        int   n    = bezf.degree();
        typedef  typename monomial_form < Point >::cpts_t cpts_t;
        cpts_t legf(sz, Point(0.0));
        double njn = 1;

        std::vector<size_t> njcns(nkcks(n,n));
        for(int k = 0; k <= n; ++k) {
            std::vector<size_t> nkjicji(nkcks(n-k, n));
            std::vector<size_t> kicis(nkcks(k,n));
            int jsign =  1;
            for(int j = 0; j <= n;++j, jsign *=- 1)
            {
                double f = sqrt(2*j + 1.0)/(n + j + 1);
                double mjk  = 0.0;
                short isign = 1;
                size_t jci  = 1;
                size_t kick = 1;
                for(int i = 0; i <= j; ++i, isign *=- 1)
                {
                    if(i!= 0) {
                        jci  *= (j - i + 1);jci /= i;
                    }
                    mjk += double(isign) * jci * kicis[i] * nkjicji[j-i];
                }
                mjk /= njcns[j];
                mjk *= f * double(jsign);
                legf[j] += mjk * b[k];
            }
        }
        double s,e;
        std::tie(s,e) = bezf.param_range();
        return legendre_form<Point>(legf,s,e);
    }
    // ../src/media/from_legendre.png

    template <class Point>
    bezier_form<Point>
    to_bezier(const legendre_form<Point>& lf)
    {
        int n  =  lf.size() - 1;
        typename bezier_form<Point>::cpts_t  cpts(n + 1, Point(0.0));

        short  ksign =  1;

        size_t ncj = 1;
        std::vector<size_t> ncjs (ncks(n));
        for(int k = 0;k <= n; ++k, ksign *=- 1)
        {
            double f   = sqrt(2*k + 1.0);

            std::vector<size_t> kcis = ncks(k);
            std::vector<size_t> nkcjs = ncks(n-k);

            for(int j = 0; j <= n; ++j)
            {
                double mjk   = 0.0;
                short  isign = 1;

                for(int i = 0;i <= j && i <= k;++i, isign *=- 1)
                {
                    if(n - k >= j  - i)
                    {
                        auto kisq = kcis[i] * kcis[i];
                        auto nkcji = nkcjs[j-i];
                        mjk += (double(isign) * kisq * nkcji);
                    }
                }
                mjk *= (ksign * f);
                mjk /= ncjs[j];
                cpts[j] +=  (mjk * lf.coeffs(k));
            }
        }
        return make_bezier_form(std::move(cpts),
                                lf.start_param(), lf.end_param());
   }
}
}
#endif // ASTI_CHANGE_BASIS

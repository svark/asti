//-*- mode:c++ -*-
#ifndef ASTI_CHANGE_BASIS
#define ASTI_CHANGE_BASIS
#include "point.hpp"
#include "bspline.hpp"
#include "bspline_queries.hpp" // for is_bezier
#include "type_utils.hpp"
#include "bezier_form.hpp"
#include "monomial_form.hpp"

namespace geom {
namespace ops  {

template<class Point>
monomial_form < Point >
to_monomial(const bspline < Point >& bezf)
{
    assert(ops::is_bezier(bezf));
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
}
}
#endif // ASTI_CHANGE_BASIS

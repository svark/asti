#pragma once
#include "tol.hpp"
namespace geom { namespace bspline_ops {
// farin cagd page 422
template <typename SplineCurve>
SplineCurve fair_by_knot_removal(SplineCurve & crv, double tol)
{
    while(crv.degree() < 3)
        crv.swap(raise_degree(crv));
    auto const &  ts =  crv.knots();
    std::vector<double> uts;
    uts.reserve(ts.size());
    std::unique_copy(ts.cbegin(),
                     ts.cend(),
                     uts.begin());
    int j = 0;
    std::vector<double> kdash_variation(uts.size());
    for(auto t:uts)
    {
        auto kdash_minus = crv.eval_derivative(3, t - tol::param_tol/ 2);
        auto kdash_plus  = crv.eval_derivative(3, t + tol::param_tol/ 2);
        kdash_variation[j] = sqlen(kdash_plus - kdash_minus);
        ++j;
    }

    auto it = std::max_element(kdash_variation.cbegin(),
                               kdash_variation.cend());
    size_t k = std::distance(kdash_variation.cbegin(), it);
    double u = uts[k];
    rmat_base_vd r(crv, u);
    size_t nu = r.locate_nu(u);
    auto & cpts = crv.control_points();
    SplineCurve::cpts_t newcpts(cpts);
    auto l =  cpts[nu - 1] * (t[nu + 1] - t[nu - 3])
        -  cpts[nu - 2] * (t[nu + 1] - t[nu]);
    l *=1.0 / (t[nu] - t[nu - 3]);
    auto r =  cpts[nu + 1] * (t[nu + 3] - t[nu - 1])
        -  cpts[nu + 2] * (t[nu] - t[nu - 1]);
    r *=1.0 / (t[nu] - t[nu - 3]);
    newcpts[nu] = (t[nu + 2] - t[nu]) * l + (t[nu] - t[nu - 2]) * r;
    newcpts[nu] *= 1 / (t[nu + 2] - t[nu - 2]);
    
    if(sqlen(newcpts[nu] - cpts[nu]) > tol * tol)
    {
        newcpts[nu] = cpts[nu] + tol *
             (newcpts[nu] - cpts[nu]) / len(newcpts[nu] - cpts[nu]);
    }
    return SplineCurve(std::move(newcpts),std::vector<double>(t), crv.degree());
}
}
}

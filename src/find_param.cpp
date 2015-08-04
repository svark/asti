#include "find_param.hpp"
#include "rational_bspline.hpp"
#include "change_basis.hpp"
#include "split_into_bezier_patches.hpp"
#include "tol.hpp"
#include <numeric>
#include "insert_knot.hpp"

namespace geom {


// root marching
std::pair<double, bool>
find_next_rootc(geom::bspline<double>& spl,
                double prev,
                double tol)
{
    auto const &t = spl.knots();
    auto const &c = spl.control_points();
    int p = spl.degree();

    size_t k = 1;
    while( k < c.size() && t[k] < prev)
        ++k;

    for(; ; )
    {
        auto const &t = spl.knots();
        auto const &c = spl.control_points();
        auto tcap = [&t, p](size_t i) -> double {
            return std::accumulate( &t[i+1], &t[i + 1] + p, 0.0 )/p;
        };

        size_t n = c.size();
        while(k < n && ( c[k-1]*c[k] > 0 ))
            ++k;

        if( k >= n)
            break;

        double root =
            tcap(k) - (c[k]/p) * (t[k+p] - t[k] )/(c[k] - c[k-1] );

        if( fabs(spl.eval(root)) < tol)
            return std::make_pair(root, true);

        ops::insert_knot(spl, root).swap(spl);
    }
    return std::make_pair(prev, false);
}


// pg 509 - 510, Hoschek
//  given a rational curve X(t) = (x(t), y(t)) where
//  x(t) =  \frac{ \sum_{i = 0}^n a_i t^i } {\sum_{i = 0}^n c_i t^i }
//  y(t) =  \frac{ \sum_{i = 0}^n b_i t^i } {\sum_{i = 0}^n c_i t^i }

//  if one writes
//  L_{ij}(x, y) =  \sum_{l \leq \text{min}(i, j), l + m = i + j + 1}
//  (b_m c_l -  c_m b_l) x +  (a_l c_m -  a_m c_l) y +  (a_m b_l - a_l
//  b_m)

// the implicit form is the Bezout resultant given the determinant of
// the matrix \left[L_{ij}\right]

// we also have \sum_j L_ij t^j =  0; for all i
// this enables us to solve for t given x, y on the curve.
// todo:sederberg suggests skipping having to get the monomial form
//
std::pair<double, bool >
find_param(point2d_t const & p,
           monomial_form < point3d_t > const & mf )
{
    std::vector<double> coeffs(mf.size());
    auto a = [&mf](int i){ return mf[i][0];};
    auto b = [&mf](int i){ return mf[i][1];};
    auto c = [&mf](int i){ return mf[i][2];};
    std::vector<double > l0j(mf.size());
    for(int j = 0;j < mf.degree(); ++j)
    {
        size_t m =  j + 1;
        l0j[j] = (b(m) * c(0) -  c(m) * b(0)) * p[0]
            + (a(0) * c(m) -  c(0) * a(m)) * p[1]
            + (a(m) * b(0) -  b(m) * a(0)) ;
    }
    bspline < double > bzf(
        ops::to_bezier(
            monomial_form < double > (l0j,
                                      mf.start_param(),
                                      mf.end_param()) ));
    return find_next_rootc(bzf, mf.start_param(), tol::resabs);
}

std::pair<double, bool >
find_param(point2d_t const & p,
           monomial_form < point2d_t > const & mf )
{
    std::vector<double> coeffs(mf.size());
    auto a = [&mf](int i){ return mf[i][0];};
    auto b = [&mf](int i){ return mf[i][1];};
    std::vector<double > l0j(mf.size());
    for(int j = 0;j < mf.degree(); ++j)
    {
        size_t m =  j + 1;
        l0j[j] = (b(m) - b(0)) * p[0]
            + (a(0) -  a(m)) * p[1]
            + (a(m) * b(0) -  b(m) * a(0));
    }
    bspline < double > bzf(
        ops::to_bezier(monomial_form < double > (l0j,
                                                 mf.start_param(),
                                                 mf.end_param()) ));
    return find_next_rootc(bzf, mf.start_param(), tol::resabs);
}

std::pair<double, bool>
find_param(point2d_t const & p,
           rational_bspline<point2d_t> const & c)
{
    auto const & patches =  ops::split_into_bezier_patches(c.spline());
    for(auto const& bp : patches)
    {
        auto const &mf = ops::to_monomial(bp);
        auto q =  find_param(p, mf);
        if(q.second && q.first <= mf.end_param())
            return q;
    }
    return std::make_pair(0, false);
}

std::pair<double, bool>
find_param(point2d_t const & p,
           bspline < point2d_t > const & c)
{
    auto const & patches =  ops::split_into_bezier_patches(c);
    for(auto const& bp : patches)
    {
        auto const &mf = ops::to_monomial(bp);
        auto q =  find_param(p, mf);
        if(q.second && q.first <= mf.end_param()) {
            if(tol::eq(len(c.eval(q.first) - p),0))
                return q;
        }
    }
    return std::make_pair(0, false);
}


}

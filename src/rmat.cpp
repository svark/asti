//-*- mode:c++ -*-
#include <cassert>
#include "geom_exception.hpp"
#include "tol.hpp"
#include "rmat.hpp"
#include <algorithm>
#include <memory>
#include <math.h>
#include "special.hpp"
#include "range.hpp"

namespace geom {

template<class KnotIter>
KnotIter
rmat_base<KnotIter>::locate(double u) const
{
    int d = deg;
    size_t n = std::distance(t, e) - d - 1;
    if(u <= t[d])
        u = (t[d] + tol::param_tol/2);
    if(u >= t[n])
        u = (t[n] - tol::param_tol/2);

    auto it = std::upper_bound(t + d, e - d, u);
    assert(it != e - d);
    return std::prev(it);
}

template<class KnotIter>
double
rmat_base<KnotIter>::der_n(size_t idx,
                           int derOrder,
                           double u) const
{
    int size = degree();
    auto nu  = locate_nu(u);
    if(idx < nu - size || idx > nu )
        return 0.0;
    std::unique_ptr<double[]> cache(new double[size+1]);
#ifdef WIN32
    std::fill_n(stdext::make_checked_array_iterator(cache.get(),size+1),
                size + 1, 0.0);
#else 
    std::fill_n(cache.get(),
		size + 1, 0.0);
#endif

    cache[idx - nu + size] = 1;
    spline_compute(nu, u, derOrder, cache.get());
    return cache[0];
}
// get the array of deg + 1 bspline basis functions at u
template<class KnotIter>
std::vector<double>
rmat_base<KnotIter>::der_coeffs_par(int derOrder, double u) const
{
    mult_rmat basis_cache;
    size_t nu = locate_nu(u);
    if(derOrder > deg)
        return basis_cache.getb();
    assert(derOrder >= 0);
    for(int j = 1;j < deg - derOrder; ++j) {
        basis_cache *= rmat_explicit<KnotIter>(t + nu, j, u);
    }
    for(int j = deg - derOrder;j <= deg; ++j) {
        basis_cache *= der_rmat_explicit<KnotIter>(t + nu, j, u);
    }
    return basis_cache.getb();
}

template<class KnotIter>
std::vector<double>
rmat_base<KnotIter>::coeffs_par(double u) const
{
    return der_coeffs_par(0, u);
}

template <class KnotIter>
Eigen::MatrixXd
rmat_base<KnotIter>::insertion_matrix(KnotIter f, KnotIter l) const
{
    size_t m = std::distance(f, l); // new knot size
    size_t n = std::distance(t, e); // old knot size
    Eigen::MatrixXd matrix(m, n);
    matrix.setZero();
    for(size_t i = 0;i < m; ++i)
    {
        mult_rmat basis_cache;
        size_t nu = locate_nu(f[i]);
        for(int j = 1;j <= deg; ++j) {
            basis_cache *= rmat_explicit<KnotIter>(t + nu, j, f[i + j]);
        }
        for(int k = 0;k < deg + 1; ++k){
            matrix(i, nu - deg + k) = basis_cache.get(k);
        }
    }
    return matrix;
}

//.  ./media/basis_comp.png
// returns the vector b = (B_{\mu - p, p}(u), \ldots, B_{\mu, p}(x)), mu =
// int such that u\in [t_\mu, t_\mu + 1) and p is the degree
// low mem foot print
template <class KnotIter>
std::vector<double>
rmat_base<KnotIter>::basis(double u)
{
    int p = degree();
    size_t nu = locate_nu(u);
    std::vector<double> b(p + 1, 0.0);
    b[p] = 1;
    for(int r = 1;r <= p; ++r)
    {
        size_t k = nu - r + 1;
        double d = (t[k + r] - t[k]);
        double lambda2 = 0;
        lambda2 = sdiv(t[k + r] - u, d);
        b[p - r] = lambda2 * b[p - r + 1];
        for(int i = p - r + 1;i < p; ++i)
        {
            ++k;
            b[i] *= (1 - lambda2);
            double d = (t[k + r] - t[k]);
            lambda2 = sdiv(t[k + r] - u, d);
            b[i] += lambda2 * b[i + 1];
        }
        b[p] *= (1 - lambda2);
    }
    return b;
}

template<class KnotIter>
size_t rmat_base<KnotIter>::start_mult() {
    double u = front();
    size_t mult = 0;
    for(auto v : util::mk_range(t, e))
    {
        if(tol::param_eq(u, v))
            ++mult;
        else
            break;
    }
    return mult;
}

template<class KnotIter>
size_t rmat_base<KnotIter>::end_mult() {
    double u =  back();
    size_t mult = 0;
    for( auto v : util::mk_rrange(t, e) )
    {
        if(tol::param_eq(u, v))
            ++mult;
        else
            break;
    }
    return mult;
}

template<class KnotIter>
size_t rmat_base<KnotIter>::mult(double u) {

    if(tol::param_eq(u , back() ))
        return end_mult();

    if(u > back() )
        return 0;

    if( tol::param_eq(u, front() ))
        return start_mult();

    if(u < front() )
        return 0;

    size_t nu = locate_nu(util::fnext(u));

    if(!tol::param_eq(u,t[nu]))
        return 0;

    size_t mult = 0;
    for( auto v : util::mk_rrange(t+nu-deg, t+nu+1) )
    {
        if(tol::param_eq(u, v))
            ++mult;
        else
            break;
    }
    return mult;
}

// tries too hard..
template <class KnotIter>
size_t rmat_base<KnotIter>::locate_nu(double u,  size_t nu_guess) const
{
    size_t size = std::distance(t, e);

    if( nu_guess + 1 < size &&
        tol::param_eq(u, t[nu_guess]))
    {
        //  we are almost there. Lets try the next strictly larger
        //  knot first.
        size_t n = nu_guess + 1;
        for(;n < size && tol::param_eq(t[n] , u); ++n)
        {
        }
        if(n < size && t[n] > u)
            return n - 1;

        else if(n == size) {       //  n is past
            size_t n = nu_guess - 1;   //  the last knot so lets try
                                    //  backwards from nu_guess
            for (;n >= 0 && tol::param_eq(t[n] , u); --n)
            {
            }
            if(n >= 0 && t[n] < u)
                return n;
        }
    }
    else if(nu_guess + 1 < size &&
            t[nu_guess] <= u && t[nu_guess + 1] > u)
    {
        return nu_guess;
    }
    return locate_nu(u);
}


template <class Point>
Point
rmat<Point>::eval_derivative(int derOrder, double u) const
{
    size_t nu = locate_nu(u);
    int size = degree();

    std::unique_ptr<point_t[]> cache(new point_t[size+1]);

    for(int j = 0; j < size + 1; ++j)
        cache[j] = control_pt(j + nu - size);

    spline_compute(nu, u, derOrder, cache.get());

    return cache[0];
}


template <class Point>
typename rmat<Point>::cpts_t
rmat<Point>::insert_knots(const std::vector<double>& taus)
{
    rmat<Point>::cpts_t newcontrol_pts;
    int size = degree();
    newcontrol_pts.reserve( taus.size() - size - 1 );

    for(size_t i = 1; i < taus.size() - size; ++i) {

        size_t nu = locate_nu(taus[i]);

        std::unique_ptr<point_t[]> cache(new point_t[size + 1]);

        for(int j = 0; j < size + 1; ++j)
            cache[j] = control_pt(j + nu - size);

        spline_compute(taus.cbegin() + i, nu, cache.get());

        newcontrol_pts.push_back(cache[0]);
    }
    return newcontrol_pts;
}
}
#include "point.hpp"
/*
Local Variables:
eval:(load-file "./scripts/temp.el")
eval:(setq methods_rmat '())
eval:(setq methods_rmat_base '())
eval:(setq point_types (list "double"
                        "point2d_t"
                        "point3d_t"
                        "point4d_t" ))
eval:(setq k_iter_types (list "const double *"
                        "std::vector<double>::iterator"
                        "std::vector<double>::const_iterator"
                        ))
eval:(instantiate-templates (file-name-sans-extension) "rmat" point_types methods_rmat (list ) )
eval:(instantiate-templates "rmat_base" "rmat_base" k_iter_types methods_rmat_base (list ) )
End:
*/
namespace geom {
// hint to the instantiation elisp that we are dealing with template classes
template <class KnotIter> struct rmat_base;
template <class Point> struct rmat;

#include "rmat_base_inst.inl"
#include "rmat_inst.inl"
}

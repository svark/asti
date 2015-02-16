//-*- mode:c++ -*-
#include <stdafx.h>
#include <cassert>
#include "geom_exception.hpp"
//#include <algorithm>
#include "tol.hpp"
#include "rmat.hpp"
#include <algorithm>
#include <memory>
#include <math.h>
#include "util.hpp"
namespace geom {


template<class KnotIter>
KnotIter
rmat_base<KnotIter>::locate(double u) const
{
    auto it = std::upper_bound(t, e, u);
    if( it == e && fabs(u - back() ) < tol::param_tol ) {
        it = std::lower_bound(t, e, back());
    }
    else if( it == t && fabs(u - front()) < tol::param_tol) {
        it = std::upper_bound(t, e, front());
    }
    if( it != e && it != t )
        return std::prev(it);

    throw geom_exception(knot_not_in_range_error);
}

template<class KnotIter>
double
rmat_base<KnotIter>::der_n(size_t idx,
                           int numDer,
                           double u) const
{
    int size = degree();
    auto nu = locate_nu(u);
    if(idx < nu - size || idx > nu )
        return 0;
    std::unique_ptr<double[]> cache(new double[size+1]);
    std::fill_n(stdext::make_checked_array_iterator(cache.get(),size+1),
                size + 1, 0.0);

    cache[idx-nu+size] = 1;
    spline_compute(nu, u, numDer, cache.get());
    return cache[0];
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

        else if(n == size ) {       //  n is past
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
rmat<Point>::eval_derivative(int numDer, double u) const
{
    size_t nu = locate_nu(u);
    int size = degree();

    std::unique_ptr<point_t[]> cache(new point_t[size+1]);

    for(size_t j = 0; j < size + 1; ++j)
        cache[j] = control_pt(j + nu - size);

    spline_compute(nu, u, numDer, cache.get());

    return cache[0];
}


template <class Point>
typename rmat<Point>::cpts_t
rmat<Point>::insert_knots(const std::vector<double>& taus)
{
    std::vector<point_t> necontrol_pts;
    int size = degree();
    necontrol_pts.reserve( taus.size() - size - 1 );

    for(size_t i = 0; i < taus.size() - size - 1; ++i) {

        size_t nu = locate_nu(taus[i]);

        std::unique_ptr<point_t[]> cache(new point_t[size + 1]);

        for(size_t j = 0; j < size + 1; ++j)
            cache[j] = control_pt(j + nu - size);

        spline_compute(taus.cbegin() + i, nu, cache.get());

        necontrol_pts.push_back(cache[0]);
    }
    return necontrol_pts;
}
}
#include "point.hpp"
/*
Local Variables:
eval:(load-file "./temp.el")
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
eval:(instantiate-templates "rmat" point_types methods_rmat (list ) )
eval:(instantiate-templates "rmat_base" k_iter_types methods_rmat_base (list ) )
End:
*/
namespace geom {
// hint to the instantiation elisp that we are dealing with template classes
template <class KnotIter> struct rmat_base;
template <class Point> struct rmat;

#include "rmat_base_inst.cpp"
#include "rmat_inst.cpp"
}

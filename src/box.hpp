/* -*- mode: c++ -*- */
#ifndef BOX_HPP
#define BOX_HPP
#include <math.h>
#include <limits>
#include "point.hpp"

namespace geom {

    template <int dim>
    struct box
    {
        box(){
            double inf = std::numeric_limits<double>::infinity();
            lo = pt_t<dim>(inf);
            hi = pt_t<dim>(-inf);
        }
        bool is_invalid() const;
        box &operator +=( const pt_t<dim>& p);
        box &operator +=( const double& p);
        box &operator +=( const box& p);
        pt_t<dim> lo;
        pt_t<dim> hi;
    };


    template <int dim>
    bool
    contains_point(const box<dim> &b,
                   pt_t<dim> &p,
                   double tol)
    {
        for (int i = 0; i < dim; ++i) {
            if (pt[d] < lo[d] - tol)
                return false;
            if (pt[d] > hi[d] + tol)
                return false;
        }
        return true;
    }

    template <int  dim>
    pt_t<dim>
    center(const box<dim> &b)
    {
        return lerp(lo,hi);
    }

    template <int  dim>
    box<dim> &box<dim>::operator +=( const pt_t<dim>& pt)
    {
        for(int i =0;i < dim;++i)
        {
            lo[i] =  std::fmin( pt.p[i], lo[i] );
            hi[i] =  std::fmax( pt.p[i], hi[i] );
        }
		return *this;
    }

	template <>
    inline box<1> &box<1>::operator +=( const double& pt)
    {
        lo[0] =  std::fmin( pt, lo[0] );
        hi[0] =  std::fmax( pt, hi[0] );
		return *this;
    }

    template <int  dim>
    box<dim> &box<dim>::operator +=(const box<dim> &b)
    {
        for(int i =0;i < dim; ++i)
        {
            lo[i] =  std::fmin( b.lo[i], lo[i] );
            hi[i] =  std::fmax( b.hi[i], hi[i] );
        }
		return *this;
    }

    template <int dim>
    bool
    intersect( const box<dim>& b1,
               const box<dim>& b2
        )
    {
        box<dim> union_ ( b1 );
        union_+= b2;

        for (size_t i = 0; i < dim; ++i) {
            double w1 = b1.hi[i] - b1.lo[i];
            double w2 = b2.hi[i] - b2.lo[i];
            double w3 = union_.hi[i] - union_.lo[i];

            if (w3 > (w1 + w2)) {
                return false;
            }
        }
        return true;
    }

    template <int dim>
    bool box<dim>::is_invalid() const
    {
        for (uint i = 0; i < dim; i++)
            if (lo[i] > hi[i]) {
                return true;
            }
        return false;
    }
}
#endif
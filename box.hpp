/* -*- mode: c++ -*- */
#ifndef BOX_HPP
#define BOX_HPP
#include <cmath>
#include <limits>
#include "point.hpp"

namespace geom {

    template <int dim>
    struct box
    {
        box(){
            double inf = std::numeric_limits<double>::infinity();
            lo = diag_point(inf);
            hi = diag_point(-inf);
        }
        bool is_invalid() const;
        box &operator +=( pt_t<dim> p);
        pt_t<dim> lo;
        pt_t<dim> hi;
    };

    template <int dim>
    struct oriented_box : public box<dim>
    {
        oriented_box(){}
        matrix3d mat;
    };

    using std::max, std::min;

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
    box &operator +=( pt_t<dim> pt)
    {
        for(int i =0;i < dim;++i)
        {
            lo[i] =  fmin( pt.p[i], lo[i] );
            hi[i] =  fmax( pt.p[i], hi[i] );
        }
    }

    template <int  dim>
    box &box::operator +=( box<dim> b)
    {
        for(int i =0;i < dim; ++i)
        {
            lo[i] =  fmin( b.lo[i], lo[i] );
            hi[i] =  fmax( b.hi[i], hi[i] );
        }
    }

    template <int dim>
    bool
    intersect( const box<dim>& b1
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

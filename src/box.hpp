/* -*- mode: c++ -*- */
#ifndef BOX_HPP
#define BOX_HPP
#include <math.h>
#include <limits>
#include "point.hpp"

namespace geom {

    template <class Point>
    struct box
    {
        box(){
            double inf = std::numeric_limits<double>::infinity();
            lo = Point(inf);
            hi = Point(-inf);
        }
        bool is_invalid() const;
        box &operator +=( const Point& p);
        box &operator +=( const box& p);
        Point lo;
        Point hi;
    };


    template <class Point>
    bool
    contains_point(const box<Point> &b,
                   Point &p,
                   double tol)
    {
        enum { dim = point_dim < Point >::dimension };
        for (int i = 0; i < dim; ++i) {
            if (pt[d] < lo[d] - tol)
                return false;
            if (pt[d] > hi[d] + tol)
                return false;
        }
        return true;
    }

    template <class Point>
    Point
    center(const box<Point> &b)
    {
        return lerp(lo,hi);
    }

    template <class Point>
    box<Point> &box<Point>::operator +=( const Point& pt)
    {
        enum { dim = point_dim < Point >::dimension };
        for(int i =0;i < dim;++i)
        {
            lo[i] =  std::fmin( pt.p[i], lo[i] );
            hi[i] =  std::fmax( pt.p[i], hi[i] );
        }
        return *this;
    }

    template <>
    inline box<double> &box<double>::operator +=( const double& pt)
    {
        lo =  std::fmin( pt, lo );
        hi =  std::fmax( pt, hi );
        return *this;
    }

    template <class Point>
    box<Point> &box<Point>::operator +=(const box<Point> &b)
    {
        enum { dim = point_dim < Point >::dimension };
        for(int i =0;i < dim; ++i)
        {
            lo[i] =  std::fmin( b.lo[i], lo[i] );
            hi[i] =  std::fmax( b.hi[i], hi[i] );
        }
        return *this;
    }

    template <class Point>
    bool
    intersect( const box<Point>& b1,
               const box<Point>& b2
        )
    {
        enum { dim = point_dim < Point >::dimension };
        box<Point> union_ ( b1 );
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

    template <class Point>
    bool box<Point>::is_invalid() const
    {
        enum { dim = point_dim < Point >::dimension };
        for (uint i = 0; i < dim; i++)
            if (lo[i] > hi[i]) {
                return true;
            }
        return false;
    }
}
#endif

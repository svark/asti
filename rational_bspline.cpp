//-*- mode:c++ -*-
#include "rational_bspline.hpp"

namespace geom
{
    //{{{ (@* "Constructors" )
    template <int dim>
    rational_bspline<dim>::rational_bspline( const cpts_t& pts,
                                            const knots_t &ks,
                                            int degree_): spl(pts,ks,degree_)
    {
    }

    template <int dim>
    rational_bspline<dim>::rational_bspline(cpts_t&& pts,
                                            knots_t &&ks,
                                            int degree_)
        : spl(std::forward(cpts),std::forward(ks),degree)
    {
    }
    template <int dim>
    rational_bspline<dim>::rational_bspline(const std::vector<point_t>& pts,
                     const std::vector<double>& weights,
                     const knots_t &ks, int degree_)
        : rational_bspline(std::move( util::interleave(pts, weights) ), ks, degree)
    {
    }

    template <int dim>
    rational_bspline<dim>::
    rational_bspline(const rational_bspline& other)
        : spl(other.spl)
    {
    }

    template <int dim>
    rational_bspline<dim>::
    rational_bspline(rational_bspline&& other)
        : spl(std::forward(other.spl))
    {
    }
    //}}}

    //{{{ (@* "Evaluators")
    template<int dim>
    vector_t
    rational_bspline<dim>::eval_derivative(double u) const
    {
        auto pt     = spl.eval(u);
        auto der_pt = spl.eval_derivative(u);

        vector_t vec;
        for(int i =0; i < dim;++i )
            vec.v[i] = ( -der_pt.v[dim] * pt.p[i] + der_pt.v[i] ) / pt.p[dim];
        return vec;
    }

    template<int dim>
    point_t
    rational_bspline<dim>::eval(double u) const
    {
        auto pt  = spl.eval(u);
        auto res = project(pt);
        return res;
    }

    template<int dim>
    template <class KnotIter>
    point_t
    rational_bspline<dim>::blossom_eval(KnotIter us) const
    {
        auto pt  = spl.blossom_eval(us);
        auto res = project(pt);
        return res;
    }
    template rational_bspline<2>::blossom_eval(double *);
    //}}}

    //{{{  (@* "split into bezier patches")
    template std::list<rational_bspline<2>> split_into_bezier_patches<>(const rational_bspline<2> &) ;
    template std::list<rational_bspline<3>> split_into_bezier_patches<>(const rational_bspline<3> &) ;
    //}}}
}

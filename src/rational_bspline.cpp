//-*- mode:c++ -*-
#include "rational_bspline.hpp"

namespace geom
{


//{{{ (@* "Constructors" )
template <class SplineType>
rational_bspline<SplineType>::rational_bspline(wcpts_t pts,
                                          knots_t ks,
                                          int degree_):
    spl(std::move(pts),std::move(ks),degree_)
{
}

template <class SplineType>
rational_bspline<SplineType>::rational_bspline(const rational_bspline& other)
    : spl(other.spl)
{
}

//}}}

//{{{ (@* "Evaluators")

//template <class SplineType>
//typename rational_bspline<SplineType>::vcpts_t
// msvc crashes with an internal error so moved this defintion to the header.
//rational_bspline<SplineType>::eval_derivatives(int numDer, double u) const
//{
//
//}


template <class SplineType>
typename rational_bspline<SplineType>::point_t
rational_bspline<SplineType>::eval(double u) const
{
    return project(spl.eval(u));
}

template <class SplineType>
template <class KnotIter>
typename rational_bspline<SplineType>::point_t
rational_bspline<SplineType>::blossom_eval(KnotIter us) const
{
    return project(spl.blossom_eval(us));
}
//}}}

}
//{{{ (@* "template instantiations")

#include "point.hpp"
#include "skip_idx_iter.hpp"
#include "bspline.hpp"
#include "periodic_bspline.hpp"
template <class SplineType> struct rational_bspline;
/*
  Local Variables:
  eval:(load-file "./scripts/temp.el")
  eval:(setq methods (list "blossom_eval") )
  eval:(setq knotIterTypes  (list "const double *"
  "std::vector<double>::const_iterator"
  "util::skip_ith_iter<std::vector<double>::const_iterator>" ) )
  eval:(setq splTypes (list "bspline<point2d_t>"
  "bspline<point3d_t>" "bspline<point4d_t>"  "periodic_bspline<point2d_t>"
  "periodic_bspline<point3d_t>" "periodic_bspline<point4d_t>"  ))
  eval:(instantiate-templates "rational_bspline"
  "rational_bspline"  splTypes  (list (cons (car methods)  knotIterTypes)) )
  End:
*/
namespace geom {
#include "rational_bspline_inst.cpp"
}
//}}}

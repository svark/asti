//-*- mode:c++ -*-
#include "rational_bspline.hpp"

namespace geom
{
//{{{ (@* "Constructors" )
template <class Point,class PTag>
rational_bspline<Point,PTag>::rational_bspline(wcpts_t pts,
                                               knots_t ks,
                                               int degree_):
    spl(std::move(pts),std::move(ks),degree_)
{
}

template <class Point,class PTag>
rational_bspline<Point,PTag>::rational_bspline(const rational_bspline& other)
    : spl(other.spl)
{
}

//}}}

//{{{ (@* "Evaluators")

//template <class Point,class PTag>
//typename rational_bspline<Point,PTag>::vcpts_t
// msvc crashes with an internal error so moved this defintion to the header.
//rational_bspline<Point,PTag>::eval_derivatives(int numDer, double u) const
//{
//
//}

template <class Point,class PTag>
typename rational_bspline<Point,PTag>::point_t
rational_bspline<Point,PTag>::eval(double u) const
{
    return project(spl.eval(u));
}

template <class Point,class PTag>
template <class KnotIter>
typename rational_bspline<Point,PTag>::point_t
rational_bspline<Point,PTag>::blossom_eval(KnotIter us) const
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
template <class Point,class PTag> struct rational_bspline;
/*
  Local Variables:
  eval:(load-file "./scripts/temp.el")
  eval:(setq methods (list "blossom_eval") )
  eval:(setq knotIterTypes  (list "const double *"
  "std::vector<double>::const_iterator"
  "util::skip_ith_iter<std::vector<double>::const_iterator>" ) )
  eval:(setq splTypes (list "point2d_t,regular_tag"
  "point3d_t,regular_tag" "double,regular_tag"  "point2d_t,periodic_tag"
  "point3d_t,periodic_tag" "double,periodic_tag"  ))
  eval:(instantiate-templates "rational_bspline"
  "rational_bspline"  splTypes  (list (cons (car methods)  knotIterTypes)) )
  End:
*/
namespace geom {
#include "rational_bspline_inst.inl"
}
//}}}

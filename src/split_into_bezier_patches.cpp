#include "tol.hpp"
#include "split_into_bezier_patches.hpp"
#include "smat.hpp"
#include <assert.h>
#include "insert_knot.hpp"
#include "Eigen/Core"
#include "rational_bspline_cons.hpp"
namespace geom {
//{{{ --(@* "split into bezier patches")

template <class SplineType>
std::list<SplineType,Eigen::aligned_allocator<SplineType> >
ops::split_into_bezier_patches(const SplineType &spl)
{
    typedef typename SplineType::knots_t knots_t;
    typedef typename SplineType::point_t point_t;

    auto &ts    = spl.knots();
    int   deg   = spl.degree();
    int   order = deg + 1;

    assert(!tol::param_eq(ts.front(), ts.back()));
    SplineType s ( clamp_end(clamp_start(spl)) );
    std::list<SplineType,Eigen::aligned_allocator<SplineType>> patches;
    knots_t newts;
    newts.reserve( s.knots().size() );
    std::unique_copy(s.knots().begin(),
                     s.knots().end(),
                     std::back_inserter(newts),
                     tol::param_eq);

    for(size_t j = 1;j < newts.size()-1; ++j)
    {
        double u = newts[j];
        patches.push_back(clamp_at_right(u, s));
        s.swap(clamp_at_left(u, s));
    }
    patches.push_back(s);
    return patches;
}

template <class SplineType>
std::list< rational_bspline < SplineType>,
           Eigen::aligned_allocator<rational_bspline < SplineType >>>
    split_into_bezier_patches(const rational_bspline < SplineType > &spl)
{
    std::list<rational_bspline < SplineType>,
              Eigen::aligned_allocator<rational_bspline < SplineType>> > res;
    auto l = split_into_bezier_patches(spl.spline());
    for(auto pc:l)
    {
        res.push_back(make_rbspline(std::move(pc)));
    }
    return res;
}


//}}}

}

/*
  Local Variables:
  eval:(load-file "./scripts/temp.el")
  eval:(setq methods (list "split_into_bezier_patches"
   ))
  eval:(setq spltypes (list "bspline<double>"
  "bspline<point2d_t>"
  "bspline<point3d_t>"
  "bspline<point4d_t>"
  "periodic_bspline<double>"
  "periodic_bspline<point2d_t>"
  "periodic_bspline<point3d_t>"
  "periodic_bspline<point4d_t>"
  "rational_bspline < bspline<point2d_t>>"
  "rational_bspline < bspline<point3d_t>>"
  "rational_bspline < bspline<point4d_t>>"
  "rational_bspline < periodic_bspline<point2d_t>>"
  "rational_bspline < periodic_bspline<point3d_t>>"
  "rational_bspline < periodic_bspline<point4d_t>>"
  ))
  eval:(instantiate-templates "split_into_bezier_patches" "ops" (list )
  (product methods spltypes ))
  End:
// dump all explicitly instantiated templates below
*/
//}}}

//{{{  instantiation
#include "bspline.hpp"
#include "periodic_bspline.hpp"
#include "point.hpp"
namespace geom {
#include "split_into_bezier_patches_inst.inl"
}
//}}}

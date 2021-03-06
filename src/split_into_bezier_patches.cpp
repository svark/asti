//-*- mode:c++ -*-
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

    assert(!tol::param_eq(spl.knots().front(), spl.knots().back()));
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
        clamp_at_left(u, s).swap(s);
    }
    patches.push_back(s);
    return patches;
}

template <class SplineType>
SplineType
ops::first_bezier_patch(const SplineType &spl)
{
    typedef typename SplineType::knots_t knots_t;

    assert(!tol::param_eq(spl.knots().front(), spl.knots().back()));
    SplineType s ( clamp_start(spl) );

    knots_t newts;
    auto ts2 = s.knots();
    assert( ts2.size() );
    newts.reserve( ts2.size() );
    double firstu= ts2.front();
    using namespace std::placeholders;
    auto biter =
        std::find_if_not(
            ts2.begin(),
            ts2.end(),
            std::bind(tol::param_eq, _1,firstu));
    assert(biter!=ts2.end());
    return clamp_at_right(*biter,s);
}

template <class SplineType>
SplineType
ops::last_bezier_patch(
    const SplineType &spl)
{
    typedef typename SplineType::knots_t knots_t;

    assert(!tol::param_eq(spl.knots().front(), spl.knots().back()));
    SplineType s(clamp_end(spl) );

    knots_t newts;
    auto ts2 = s.knots();
    assert( ts2.size() );
    newts.reserve( ts2.size() );
    double lastu= ts2.back();
    using namespace std::placeholders;
    auto biter =
        std::find_if_not(
            ts2.rbegin(),
            ts2.rend(),
            std::bind(tol::param_eq, _1,lastu));
    assert(biter!=ts2.rend());
    return clamp_at_left(*biter,s);
}

//}}}

}

/*
  Local Variables:
  eval:(load-file "./scripts/temp.el")
  eval:(setq methods (list "split_into_bezier_patches"  "first_bezier_patch" "last_bezier_patch"
  ))
  eval:(setq spltypes (list "bspline<double>"
  "bspline<point2d_t>"
  "bspline<point3d_t>"
  "bspline<point4d_t>"
  "rational_bspline<double>"
  "rational_bspline<point2d_t>"
  "rational_bspline<point3d_t>"
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

//-*- mode:c++ -*-
#ifndef ASTI_MODIFY_BSPLINE
#define ASTI_MODIFY_BSPLINE
#include "bspline_fwd.hpp"
namespace geom { namespace ops {

// Fn should be a function of type
//    Fn(SplineType::cpts_t::iterator b,e,
//       SplineType::knots_t::iterator tsb,tse)
//  This allows tweaking the values of individual control points and/or knots
template <class SplineType, class Fn>
SplineType&
modify_bspline(SplineType& spl, Fn f)
{
    SplineType(std::move(spl),f).swap(spl);
    return spl;
}

}
}
#endif // ASTI_MODIFY_BSPLINE

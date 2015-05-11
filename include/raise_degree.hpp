#ifndef ASTI_RAISE_DEGREE
#define ASTI_RAISE_DEGREE
#include "bspline_fwd.hpp"
#include "rational_bspline_cons.hpp"
namespace geom {
namespace ops
{

template <class SplineType>
extern SplineType
raise_degree(const SplineType&crv);

template<class SplineCurve>
bool
match_degrees(SplineCurve & spl1,
              SplineCurve & spl2)
{
    int p1 = spl1.degree();
    int p2 = spl2.degree();

    for(;p2 < p1;++p2)
      raise_degree(spl2).swap(spl2);

    for(;p1 < p2;++p1)
      raise_degree(spl1).swap(spl1);
    return true;
}
}
}
#endif // ASTI_RAISE_DEGREE

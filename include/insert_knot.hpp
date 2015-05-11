#ifndef ASTI_INSERT_KNOT
#define ASTI_INSERT_KNOT
#include <vector>
#include "rational_bspline_cons.hpp"
namespace geom{
namespace ops {

template <class SplineType>
extern  SplineType
insert_knots(const SplineType& crv,
             const std::vector<double>& refined_knots);

template <class SplineType>
extern  SplineType
insert_knot(const SplineType& crv, double u);


template <class SplineCurve>
bool
match_knots(SplineCurve & spl1, SplineCurve & spl2)
{
    assert(spl1.degree() == spl2.degree());
    assert(is_regular(spl1) && is_regular(spl2));
    auto & rspl1 = reparametrize(spl1); // get a 0, 1 parametrization 
    auto & rspl2 = reparametrize(spl2);
    std::vector<double> t;
    t.reserve(rspl1.knots().size()
              +  rspl2.knots().size());

    std::merge(rspl1.knots().cbegin(), rspl1.knots().cend(),
               rspl2.knots().cbegin(), std::back_inserter(t));
	auto & ispl1 = insert_knots(rspl1, t);
	auto & ispl2 = insert_knots(rspl2, t);
    spl1.swap(ispl1);
    spl2.swap(ispl2);
    return true;
}

}
}
#endif // ASTI_INSERT_KNOT

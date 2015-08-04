#ifndef ASTI_FIND_PARAM
#define ASTI_FIND_PARAM
#include "bspline_fwd.hpp"
#include "point_fwd.hpp"
#include "spline_traits.hpp"
#include <utility>
namespace geom {

extern
std::pair<double, bool>
find_next_rootc(geom::bspline<double>& spl,
                double prev,
                double tol);

extern std::pair<double, bool>
find_param(point2d_t const & p,
           geom::rational_bspline < point2d_t, regular_tag> const & c );

extern std::pair<double, bool>
find_param(point2d_t const & p,
           geom::bspline < point2d_t > const & c );
}
#endif // ASTI_FIND_PARAM

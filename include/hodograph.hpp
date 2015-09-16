#ifndef ASTI_HODOGRAPH
#define ASTI_HODOGRAPH
#include "bspline_fwd.hpp"

namespace geom{ namespace qry {

template <class Point>
extern bspline<Point> hodograph(const bspline<Point>& spl,int num_der = 1);

}
}
#endif // ASTI_HODOGRAPH

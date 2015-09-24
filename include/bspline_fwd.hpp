#ifndef ASTI_BSPLINE_FWD
#define ASTI_BSPLINE_FWD
namespace geom {
template <class Point> class bspline;
template <class Point> class periodic_bspline;
template <class Point,class PTag> class rational_bspline ;
template <class Point>
bool check_invariants(const bspline<Point>& spl);
template <class Point,class PTag>
bool check_invariants(const rational_bspline<Point,PTag>& spl);
template <class Point>
bool check_invariants(const periodic_bspline<Point>& spl);

}
#endif // ASTI_BSPLINE_FWD

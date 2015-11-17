#ifndef ASTI_BSPLINE_QUERIES
#define ASTI_BSPLINE_QUERIES
#include "bspline_fwd.hpp"
#include "tol.hpp"
#include "spline_traits.hpp"
#include "point_dim.hpp"
#include <limits>
namespace geom {
namespace qry {

template <class SplineType>
extern bool is_bezier(const SplineType& spl);

template <class Point,class PTag>
bool is_bezier(const rational_bspline < Point, PTag > & spl)
{
    return is_bezier(spl.spline());
}

template <class SplineType>
extern bool is_periodic(const SplineType & spl);

template <class SplineType>
extern bool is_clamped(const SplineType & spl);

template <class Point,class PTag>
bool is_periodic(const rational_bspline < Point, PTag >& spl)
{
    return is_periodic(spl.spline());
}

template <class Point,class PTag>
bool is_clamped(const rational_bspline < Point,PTag >& spl)
{
    return is_clamped(spl.spline());
}

template <class SplineType>
bool is_open(const  SplineType & spl)
{
    return !is_clamped(spl);
}

template <class Point>
double weight(const Point&p, rational_tag )
{
    return p[point_dim<Point>::dimension-1];
}

template <class Point>
double weight(const Point&p, polynomial_tag )
{
    return 1;
}

template <class Point>
auto lift_dim(const Point& p1, polynomial_tag) -> RAWTYPE(higher_dim(p1))
{
    typedef RAWTYPE(higher_dim(p1)) PointW ;
    return PointW(p1);
}

template <class Point>
auto lift_dim(const Point& p1, rational_tag) -> RAWTYPE(higher_dim(p1))
{
    typedef RAWTYPE(higher_dim(p1)) PointW ;
    PointW p2(p1);
    enum {dim = point_dim<PointW>::dimension };
    std::swap(p2[dim],p2[dim-1]); // move the weight to the last coordinate
    return p2;
}

template <class SplineType>
extern double greville(const SplineType& spl,int i);

extern point3d_t
auto_lift_dim3(const point3d_t& p1, polynomial_tag, polynomial_tag);

extern point3d_t auto_lift_dim3(const point2d_t& p1, polynomial_tag, polynomial_tag);

extern point4d_t
auto_lift_dim3(const point2d_t& p1, polynomial_tag, rational_tag);

extern point4d_t
auto_lift_dim3(const point3d_t& p1, polynomial_tag, rational_tag);

extern point4d_t
auto_lift_dim3(const point3d_t& p1, rational_tag, rational_tag);

extern point4d_t
auto_lift_dim3(const point4d_t& p1, rational_tag, rational_tag);


extern point3d_t
auto_lift_dim3(const point4d_t& p1, rational_tag, polynomial_tag);


extern point3d_t
auto_lift_dim3(const point3d_t& p1, rational_tag, polynomial_tag);


namespace detail {
template <class Spl>
auto get_spline(const Spl& s, rational_tag, periodic_tag)
{
    return s.spline().spline();
}

template <class Spl>
auto get_spline(const Spl& s, polynomial_tag, periodic_tag)
{
    return s.spline();
}

template <class Spl>
auto get_spline(const Spl& s, rational_tag, regular_tag)
{
    return s.spline();
}

template <class Spl>
//regular bspline
const Spl& get_spline(const Spl& s, polynomial_tag, regular_tag)
{
    return s;
}
}

template <class Spl>
auto get_spline(const Spl& s)
{
    return detail::get_spline(
        s,
        typename spline_traits<Spl>::rtag(),
        typename spline_traits<Spl>::ptag()
        );
}

template <class Crv>
size_t num_knots(const Crv& spl)
{
    return get_spline(spl).knots().size();
}

template <class Crv>
size_t num_cpts(const Crv& spl)
{
    return get_spline(spl).control_points().size();
}

template <class SplineType>
extern double curvature(const SplineType& spl, double u);

template <class SplineType>
extern double torsion(const SplineType & spl, double u);

template <class Fn, class SplineType>
auto checked_op(Fn f, const SplineType& spl, double u) -> decltype(f(spl,u))
{
    auto v = f(spl, u);
    if(std::isnan(coord(v,X)))
    {
        auto bef = f(spl, u - tol::param_tol);
        auto aft = f(spl, u + tol::param_tol);
        return (bef+aft)/2.0;
    }
    return v;
}

template <class SplineType>
double curvature_approx(const SplineType& spl, double u)
{
    return checked_op( & curvature<SplineType>, spl, u);
}

template <class SplineType>
double torsion_approx(const SplineType & spl, double u)
{
    return checked_op( & torsion<SplineType>, spl, u);
}

template <class SplineType>
auto
derivative_approx(const SplineType& spl, double u, int der)
{
    auto op = [der](const SplineType& spl, double u)
        { return spl.eval_derivative(der,u); };

    return checked_op(op, spl, u);
}



template <class SplineType>
std::pair<double,double>
param_range(const SplineType&spl)
{
    return spl.param_range();
}

template <class SplineType>
double
start_param(const SplineType&spl)
{
    return spl.param_range().first;
}
template <class SplineType>
double
end_param(const SplineType&spl)
{
    return spl.param_range().second;
}

}}

#endif

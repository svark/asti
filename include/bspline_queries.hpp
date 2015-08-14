#ifndef ASTI_BSPLINE_QUERIES
#define ASTI_BSPLINE_QUERIES
#include "bspline_fwd.hpp"
#include "tol.hpp"
#include "geom_exception.hpp"
#include "spline_traits.hpp"

namespace geom {
namespace ops {

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
double weight(const Point&p, rational_tag,
              std::integral_constant<int,3> )
{
    return 1;
}

template <class Point>
double weight(const Point&p, rational_tag,
              std::integral_constant<int,2> )
{
    return 1;
}

inline double weight(const point4d_t&p, rational_tag,
              std::integral_constant<int,3> spaceDim)
{
    return coord(p,spaceDim);
}

inline double weight(const point3d_t&p, rational_tag,
              std::integral_constant<int,2> spaceDim)
{
    return coord(p,spaceDim);
}

template <class Point>
double weight(const Point&p, polynomial_tag, std::integral_constant<int,2>)
{
    return 1;
}

template <class Point>
double weight(const Point&p, polynomial_tag, std::integral_constant<int,3>)
{
    return 1;
}
template <class Point1,class Point2>
void
lift_dim(const Point1& p1, polynomial_tag,   Point2& p2)
{
	p2 = Point2(p1);
}

template <class Point1>
void
lift_dim(const Point1& p1, rational_tag,   Point1& p2)
{
	p2 = Point1(p1);
}

template <int dim>
void 
lift_dim(const pt_t<dim>& p1, rational_tag, pt_t<dim+1>& p2)
{
	p2 = pt_t<dim+1>(pt_t<dim-1>(p1), 0.0);
	p2[dim] = p1[dim-1];
}

namespace detail {
template <class Spl>
auto get_spline(const Spl& s, std::true_type)
    -> decltype(s.spline())
{
    return s.spline();
}

template <class Spl>
const Spl& get_spline(const Spl& s, std::false_type)
{
    return s;
}
}

template <class Spl>
auto get_spline(const Spl& s) -> decltype(
    detail::get_spline(s,
                       std::integral_constant<
                       bool,
                       is_rational_type<Spl>::value
                       || is_periodic_type<Spl>::value >()) )
{
    return detail::get_spline(
        s,std::integral_constant<
        bool,is_rational_type<Spl>::value
        || is_periodic_type<Spl>::value >());
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
double checked_op(Fn f, const SplineType& spl, double u)
{
    try {
        f(spl, u);
    }catch(geom_exception e)
    {
        if(e.code() != knot_not_in_range_error_der)
            throw;

        return (f(spl, u - tol::param_tol) +
                f(spl, u + tol::param_tol)) / 2;
    }
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

}}

#endif

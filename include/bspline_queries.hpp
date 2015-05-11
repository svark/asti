#include "bspline_fwd.hpp"
#include "tol.hpp"
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
extern bool is_regular(const SplineType & spl);

template <class Point,class PTag>
bool is_periodic(const rational_bspline < Point, PTag >& spl)
{
    return is_periodic(spl.spline());
}

template <class Point,class PTag>
bool is_regular(const rational_bspline < Point,PTag >& spl)
{
    return is_regular(spl.spline());
}

template <class SplineType>
bool is_open(const  SplineType & spl)
{
    return !is_regular(spl);
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

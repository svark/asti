#ifndef ASTI_IMPLICIT
#define ASTI_IMPLICIT
#include <vector>
#include <memory>
#include "point_fwd.hpp"
#include "bspline_fwd.hpp"
#include "spline_traits.hpp"

namespace geom {
struct implicitCurveFormBase
{
    implicitCurveFormBase(std::vector<double> coefficients_,
                          int qdeg_):coefficients(std::move(coefficients_)),
        qdeg(qdeg_)
    {
    }
	double eval3d(const point3d_t &p3d) const;
    virtual double eval(const point3d_t& p) const = 0 ;
	double eval(const point2d_t& p) const;
protected:
    int qdeg;
    std::vector<double> coefficients;
};


std::unique_ptr < implicitCurveFormBase >
implicitize(const rational_bspline<point2d_t, regular_tag>& spl, int qdeg);

std::unique_ptr < implicitCurveFormBase >
implicitize(const bspline<point2d_t>& spl, int qdeg);
}
#endif // ASTI_IMPLICIT

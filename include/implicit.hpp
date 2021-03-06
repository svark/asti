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
    virtual ~implicitCurveFormBase(){}
protected:
    std::vector<double> coefficients;
    int qdeg;
};


template <class PTag>
extern std::vector<std::unique_ptr < implicitCurveFormBase> >
implicitize(const rational_bspline<point2d_t, PTag>& spl, int qdeg);

extern std::vector<std::unique_ptr < implicitCurveFormBase> >
implicitize(const bspline<point2d_t>& spl, int qdeg);

extern std::unique_ptr < implicitCurveFormBase >
implicitize(const std::function<point3d_t(double )>&, int qdeg, int sdeg);

}
#endif // ASTI_IMPLICIT

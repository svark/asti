#ifndef ASTI_TOL_HPP
#define ASTI_TOL_HPP
#include <math.h> //fabs
#include "point_fwd.hpp"
namespace tol
{
static const double param_tol = 1e-10;
static const double resabs    = 1e-8;
static const double sqresabs  = 1e-16;

inline bool eq(double v,double w, double tol=resabs) { return fabs(v-w)<tol; }
inline bool neq(double v,double w, double tol=resabs) { return fabs(v-w)>=tol; }

inline bool eqres(double v,double w) { return fabs(v-w)<resabs; }
inline bool neqres(double v,double w) { return fabs(v-w)>=resabs; }

inline bool param_eq(double v,double w) { return fabs(v-w)<param_tol; }
inline bool not_small(double v, double tol = resabs) { return fabs(v) >= tol; }
inline bool small(double v, double tol = resabs) { return fabs(v) < tol; }
template <class Point>
bool pt_eq(const Point& p1, const Point& p2) { return geom::len(p1-p2) < resabs; }
}
#endif //TOL_HPP

#ifndef TOL_HPP
#define TOL_HPP
#include <math.h>
namespace tol
{
  static double param_tol= 1e-10;
  static double resabs   = 1e-8;

  inline bool eq(double v,double w, double tol=resabs) { return fabs(v-w)<tol; }
  inline bool param_eq(double v,double w) { return fabs(v-w)<param_tol; }

}
#endif //TOL_HPP

namespace tol
{
  static double param_tol= 1e-10;
  static double resabs   = 1e-8;
  double eq(double v,double w, double tol) {
      return fabs(v-w) < tol;
  }
}

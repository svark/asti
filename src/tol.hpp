namespace tol
{
  static double param_tol= 1e-10;
  static double resabs   = 1e-8;

  inline double eq(double v,double w, double tol) { return fabs(v-w)<tol; }
  inline double param_eq(double v,double w) { return fabs(v-w)<param_tol; }

}

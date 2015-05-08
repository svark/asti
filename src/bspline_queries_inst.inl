//-*-mode:c++-*-
//Generated on: Fri Apr 17 12:36:45 2015. Do not edit
//________________________________________________________
// method:is_periodic
template bool ops::is_periodic(const bspline<double> & crv) ;
template bool ops::is_periodic(const bspline<point2d_t> & crv) ;
template bool ops::is_periodic(const bspline<point3d_t> & crv) ;
template bool ops::is_periodic(const bspline<point4d_t> & crv) ;
//________________________________________________________
// method:is_bezier
template bool  ops::is_bezier(const bspline<double>& spl) ;
template bool  ops::is_bezier(const bspline<point2d_t>& spl) ;
template bool  ops::is_bezier(const bspline<point3d_t>& spl) ;
template bool  ops::is_bezier(const bspline<point4d_t>& spl) ;
//________________________________________________________
// method:curvature
template double  ops::curvature(const bspline<double> & c, double u) ;
template double  ops::curvature(const bspline<point2d_t> & c, double u) ;
template double  ops::curvature(const bspline<point3d_t> & c, double u) ;
template double  ops::curvature(const bspline<point4d_t> & c, double u) ;
//________________________________________________________
// method:torsion
template double  ops::torsion(const bspline<double> & spl, double u) ;
template double  ops::torsion(const bspline<point2d_t> & spl, double u) ;
template double  ops::torsion(const bspline<point3d_t> & spl, double u) ;
template double  ops::torsion(const bspline<point4d_t> & spl, double u) ;

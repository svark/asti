//-*-mode:c++-*-
//Generated on: Fri Aug 14 17:57:51 2015. Do not edit
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
// method:is_clamped
template bool  ops::is_clamped(const bspline<double> & c) ;
template bool  ops::is_clamped(const bspline<point2d_t> & c) ;
template bool  ops::is_clamped(const bspline<point3d_t> & c) ;
template bool  ops::is_clamped(const bspline<point4d_t> & c) ;
//________________________________________________________
// method:curvature
template double  ops::curvature(const bspline<double> & spl, double u) ;
template double  ops::curvature(const bspline<point2d_t> & spl, double u) ;
template double  ops::curvature(const bspline<point3d_t> & spl, double u) ;
template double  ops::curvature(const bspline<point4d_t> & spl, double u) ;
//________________________________________________________
// method:torsion
template double  ops::torsion(const bspline<double> & spl, double u) ;
template double  ops::torsion(const bspline<point2d_t> & spl, double u) ;
template double  ops::torsion(const bspline<point3d_t> & spl, double u) ;
template double  ops::torsion(const bspline<point4d_t> & spl, double u) ;

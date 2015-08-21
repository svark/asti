//-*-mode:c++-*-
//Generated on: Fri Aug 21 11:20:02 2015. Do not edit
//________________________________________________________
// method:is_periodic
template bool qry::is_periodic(const bspline<double> & crv) ;
template bool qry::is_periodic(const bspline<point2d_t> & crv) ;
template bool qry::is_periodic(const bspline<point3d_t> & crv) ;
template bool qry::is_periodic(const bspline<point4d_t> & crv) ;
//________________________________________________________
// method:is_bezier
template bool  qry::is_bezier(const bspline<double>& spl) ;
template bool  qry::is_bezier(const bspline<point2d_t>& spl) ;
template bool  qry::is_bezier(const bspline<point3d_t>& spl) ;
template bool  qry::is_bezier(const bspline<point4d_t>& spl) ;
//________________________________________________________
// method:is_clamped
template bool  qry::is_clamped(const bspline<double> & c) ;
template bool  qry::is_clamped(const bspline<point2d_t> & c) ;
template bool  qry::is_clamped(const bspline<point3d_t> & c) ;
template bool  qry::is_clamped(const bspline<point4d_t> & c) ;
//________________________________________________________
// method:curvature
template double  qry::curvature(const bspline<double> & spl, double u) ;
template double  qry::curvature(const bspline<point2d_t> & spl, double u) ;
template double  qry::curvature(const bspline<point3d_t> & spl, double u) ;
template double  qry::curvature(const bspline<point4d_t> & spl, double u) ;
//________________________________________________________
// method:torsion
template double  qry::torsion(const bspline<double> & spl, double u) ;
template double  qry::torsion(const bspline<point2d_t> & spl, double u) ;
template double  qry::torsion(const bspline<point3d_t> & spl, double u) ;
template double  qry::torsion(const bspline<point4d_t> & spl, double u) ;

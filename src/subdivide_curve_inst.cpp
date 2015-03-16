//-*-mode:c++-*-
//Generated on: Wed Feb 11 14:35:13 2015. Do not edit
//________________________________________________________
// method:split_open_curve
template std::pair<bspline<double>, bspline<double>> bspline_ops::split_open_curve(const bspline<double>& spl, double u) ;
template std::pair<bspline<point2d_t>, bspline<point2d_t>> bspline_ops::split_open_curve(const bspline<point2d_t>& spl, double u) ;
template std::pair<bspline<point3d_t>, bspline<point3d_t>> bspline_ops::split_open_curve(const bspline<point3d_t>& spl, double u) ;
template std::pair<bspline<point4d_t>, bspline<point4d_t>> bspline_ops::split_open_curve(const bspline<point4d_t>& spl, double u) ;
//________________________________________________________
// method:split_periodic_curve
template bspline<double> bspline_ops::split_periodic_curve(const periodic_bspline<double>& pspl, double u) ;
template bspline<point2d_t> bspline_ops::split_periodic_curve(const periodic_bspline<point2d_t>& pspl, double u) ;
template bspline<point3d_t> bspline_ops::split_periodic_curve(const periodic_bspline<point3d_t>& pspl, double u) ;
template bspline<point4d_t> bspline_ops::split_periodic_curve(const periodic_bspline<point4d_t>& pspl, double u) ;

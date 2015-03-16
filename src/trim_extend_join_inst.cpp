//-*-mode:c++-*-
//Generated on: Thu Feb 12 10:12:50 2015. Do not edit
//________________________________________________________
// method:extract_regular_curve
template bspline<double> bspline_ops::extract_regular_curve(const bspline<double> &spl) ;
template bspline<point2d_t> bspline_ops::extract_regular_curve(const bspline<point2d_t> &spl) ;
template bspline<point3d_t> bspline_ops::extract_regular_curve(const bspline<point3d_t> &spl) ;
template bspline<point4d_t> bspline_ops::extract_regular_curve(const bspline<point4d_t> &spl) ;
//________________________________________________________
// method:trim_curve
template bspline<double> bspline_ops::trim_curve(const bspline<double> &spl, double a, double b) ;
template bspline<point2d_t> bspline_ops::trim_curve(const bspline<point2d_t> &spl, double a, double b) ;
template bspline<point3d_t> bspline_ops::trim_curve(const bspline<point3d_t> &spl, double a, double b) ;
template bspline<point4d_t> bspline_ops::trim_curve(const bspline<point4d_t> &spl, double a, double b) ;
//________________________________________________________
// method:extend_curve_start
template bspline<double> bspline_ops::extend_curve_start(const bspline<double> & spl,
                                double delta) ;
template bspline<point2d_t> bspline_ops::extend_curve_start(const bspline<point2d_t> & spl,
                                double delta) ;
template bspline<point3d_t> bspline_ops::extend_curve_start(const bspline<point3d_t> & spl,
                                double delta) ;
template bspline<point4d_t> bspline_ops::extend_curve_start(const bspline<point4d_t> & spl,
                                double delta) ;
//________________________________________________________
// method:extend_curve_end
template bspline<double> bspline_ops::extend_curve_end(const bspline<double> & spl,
                              double delta) ;
template bspline<point2d_t> bspline_ops::extend_curve_end(const bspline<point2d_t> & spl,
                              double delta) ;
template bspline<point3d_t> bspline_ops::extend_curve_end(const bspline<point3d_t> & spl,
                              double delta) ;
template bspline<point4d_t> bspline_ops::extend_curve_end(const bspline<point4d_t> & spl,
                              double delta) ;
//________________________________________________________
// method:join_starts
template bspline<double> bspline_ops::join_starts(const bspline<double>& spl1,
                         const bspline<double>& spl2,
                         int join_cont) ;
template bspline<point2d_t> bspline_ops::join_starts(const bspline<point2d_t>& spl1,
                         const bspline<point2d_t>& spl2,
                         int join_cont) ;
template bspline<point3d_t> bspline_ops::join_starts(const bspline<point3d_t>& spl1,
                         const bspline<point3d_t>& spl2,
                         int join_cont) ;
template bspline<point4d_t> bspline_ops::join_starts(const bspline<point4d_t>& spl1,
                         const bspline<point4d_t>& spl2,
                         int join_cont) ;

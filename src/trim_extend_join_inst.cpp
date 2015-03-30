//-*-mode:c++-*-
//Generated on: Fri Mar 27 20:18:42 2015. Do not edit
//________________________________________________________
// method:extract_regular_curve
template bspline<double> ops::extract_regular_curve(const bspline<double> &spl) ;
template bspline<point2d_t> ops::extract_regular_curve(const bspline<point2d_t> &spl) ;
template bspline<point3d_t> ops::extract_regular_curve(const bspline<point3d_t> &spl) ;
template bspline<point4d_t> ops::extract_regular_curve(const bspline<point4d_t> &spl) ;
template rational_bspline < bspline < point2d_t >> ops::extract_regular_curve(const rational_bspline < bspline < point2d_t >> &spl) ;
template rational_bspline < bspline < point3d_t >> ops::extract_regular_curve(const rational_bspline < bspline < point3d_t >> &spl) ;
template rational_bspline < bspline < point4d_t >> ops::extract_regular_curve(const rational_bspline < bspline < point4d_t >> &spl) ;
//________________________________________________________
// method:trim_curve
template bspline<double> ops::trim_curve(const bspline<double> &spl, double a, double b) ;
template bspline<point2d_t> ops::trim_curve(const bspline<point2d_t> &spl, double a, double b) ;
template bspline<point3d_t> ops::trim_curve(const bspline<point3d_t> &spl, double a, double b) ;
template bspline<point4d_t> ops::trim_curve(const bspline<point4d_t> &spl, double a, double b) ;
template rational_bspline < bspline < point2d_t >> ops::trim_curve(const rational_bspline < bspline < point2d_t >> &spl, double a, double b) ;
template rational_bspline < bspline < point3d_t >> ops::trim_curve(const rational_bspline < bspline < point3d_t >> &spl, double a, double b) ;
template rational_bspline < bspline < point4d_t >> ops::trim_curve(const rational_bspline < bspline < point4d_t >> &spl, double a, double b) ;
//________________________________________________________
// method:extend_curve_start
template bspline<double> ops::extend_curve_start(const bspline<double> & spl, double delta) ;
template bspline<point2d_t> ops::extend_curve_start(const bspline<point2d_t> & spl, double delta) ;
template bspline<point3d_t> ops::extend_curve_start(const bspline<point3d_t> & spl, double delta) ;
template bspline<point4d_t> ops::extend_curve_start(const bspline<point4d_t> & spl, double delta) ;
template rational_bspline < bspline < point2d_t >> ops::extend_curve_start(const rational_bspline < bspline < point2d_t >> & spl, double delta) ;
template rational_bspline < bspline < point3d_t >> ops::extend_curve_start(const rational_bspline < bspline < point3d_t >> & spl, double delta) ;
template rational_bspline < bspline < point4d_t >> ops::extend_curve_start(const rational_bspline < bspline < point4d_t >> & spl, double delta) ;
//________________________________________________________
// method:extend_curve_end
template bspline<double> ops::extend_curve_end(const bspline<double> & spl, double delta) ;
template bspline<point2d_t> ops::extend_curve_end(const bspline<point2d_t> & spl, double delta) ;
template bspline<point3d_t> ops::extend_curve_end(const bspline<point3d_t> & spl, double delta) ;
template bspline<point4d_t> ops::extend_curve_end(const bspline<point4d_t> & spl, double delta) ;
template rational_bspline < bspline < point2d_t >> ops::extend_curve_end(const rational_bspline < bspline < point2d_t >> & spl, double delta) ;
template rational_bspline < bspline < point3d_t >> ops::extend_curve_end(const rational_bspline < bspline < point3d_t >> & spl, double delta) ;
template rational_bspline < bspline < point4d_t >> ops::extend_curve_end(const rational_bspline < bspline < point4d_t >> & spl, double delta) ;
//________________________________________________________
// method:join_starts
template bspline<double> ops::join_starts(const bspline<double>& spl1,
                         const bspline<double>& spl2,
                         int join_cont) ;
template bspline<point2d_t> ops::join_starts(const bspline<point2d_t>& spl1,
                         const bspline<point2d_t>& spl2,
                         int join_cont) ;
template bspline<point3d_t> ops::join_starts(const bspline<point3d_t>& spl1,
                         const bspline<point3d_t>& spl2,
                         int join_cont) ;
template bspline<point4d_t> ops::join_starts(const bspline<point4d_t>& spl1,
                         const bspline<point4d_t>& spl2,
                         int join_cont) ;
template rational_bspline < bspline < point2d_t >> ops::join_starts(const rational_bspline < bspline < point2d_t >>& spl1,
                         const rational_bspline < bspline < point2d_t >>& spl2,
                         int join_cont) ;
template rational_bspline < bspline < point3d_t >> ops::join_starts(const rational_bspline < bspline < point3d_t >>& spl1,
                         const rational_bspline < bspline < point3d_t >>& spl2,
                         int join_cont) ;
template rational_bspline < bspline < point4d_t >> ops::join_starts(const rational_bspline < bspline < point4d_t >>& spl1,
                         const rational_bspline < bspline < point4d_t >>& spl2,
                         int join_cont) ;
//________________________________________________________
// method:extend_curve_end_to_pt
template bspline<double> ops::extend_curve_end_to_pt(const bspline<double> & spl,
                                     bspline<double>::point_t const & target) ;
template bspline<point2d_t> ops::extend_curve_end_to_pt(const bspline<point2d_t> & spl,
                                     bspline<point2d_t>::point_t const & target) ;
template bspline<point3d_t> ops::extend_curve_end_to_pt(const bspline<point3d_t> & spl,
                                     bspline<point3d_t>::point_t const & target) ;
template bspline<point4d_t> ops::extend_curve_end_to_pt(const bspline<point4d_t> & spl,
                                     bspline<point4d_t>::point_t const & target) ;
template rational_bspline < bspline < point2d_t >> ops::extend_curve_end_to_pt(const rational_bspline < bspline < point2d_t >> & spl,
                                     rational_bspline < bspline < point2d_t >>::point_t const & target) ;
template rational_bspline < bspline < point3d_t >> ops::extend_curve_end_to_pt(const rational_bspline < bspline < point3d_t >> & spl,
                                     rational_bspline < bspline < point3d_t >>::point_t const & target) ;
template rational_bspline < bspline < point4d_t >> ops::extend_curve_end_to_pt(const rational_bspline < bspline < point4d_t >> & spl,
                                     rational_bspline < bspline < point4d_t >>::point_t const & target) ;

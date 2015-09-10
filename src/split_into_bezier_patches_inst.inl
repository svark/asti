//-*-mode:c++-*-
//Generated on: Tue Sep  8 20:08:42 2015. Do not edit
//________________________________________________________
// method:split_into_bezier_patches
template std::list<bspline<double>,Eigen::aligned_allocator<bspline<double>> > ops::split_into_bezier_patches(const bspline<double> &spl) ;
template std::list<bspline<point2d_t>,Eigen::aligned_allocator<bspline<point2d_t>> > ops::split_into_bezier_patches(const bspline<point2d_t> &spl) ;
template std::list<bspline<point3d_t>,Eigen::aligned_allocator<bspline<point3d_t>> > ops::split_into_bezier_patches(const bspline<point3d_t> &spl) ;
template std::list<bspline<point4d_t>,Eigen::aligned_allocator<bspline<point4d_t>> > ops::split_into_bezier_patches(const bspline<point4d_t> &spl) ;
//________________________________________________________
// method:first_bezier_patch
template bspline<double> ops::first_bezier_patch(const bspline<double> &spl) ;
template bspline<point2d_t> ops::first_bezier_patch(const bspline<point2d_t> &spl) ;
template bspline<point3d_t> ops::first_bezier_patch(const bspline<point3d_t> &spl) ;
template bspline<point4d_t> ops::first_bezier_patch(const bspline<point4d_t> &spl) ;
//________________________________________________________
// method:last_bezier_patch
template bspline<double> ops::last_bezier_patch(const bspline<double> &spl) ;
template bspline<point2d_t> ops::last_bezier_patch(const bspline<point2d_t> &spl) ;
template bspline<point3d_t> ops::last_bezier_patch(const bspline<point3d_t> &spl) ;
template bspline<point4d_t> ops::last_bezier_patch(const bspline<point4d_t> &spl) ;

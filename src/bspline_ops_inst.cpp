//-*-mode:c++-*-
//Generated on: Mon Feb  2 18:11:15 2015. Do not edit
//________________________________________________________
// method:insert_knots
template bspline<double> bspline_ops::insert_knots(const bspline<double>& crv,
                          const std::vector<double>& new_knots);
template bspline<point2d_t> bspline_ops::insert_knots(const bspline<point2d_t>& crv,
                          const std::vector<double>& new_knots);
template bspline<point3d_t> bspline_ops::insert_knots(const bspline<point3d_t>& crv,
                          const std::vector<double>& new_knots);
template bspline<point4d_t> bspline_ops::insert_knots(const bspline<point4d_t>& crv,
                          const std::vector<double>& new_knots);
//________________________________________________________
// method:reparametrize
template bspline<double>  bspline_ops::reparametrize(const bspline<double>& spl,
                                             double t1, double t2);
template bspline<point2d_t>  bspline_ops::reparametrize(const bspline<point2d_t>& spl,
                                             double t1, double t2);
template bspline<point3d_t>  bspline_ops::reparametrize(const bspline<point3d_t>& spl,
                                             double t1, double t2);
template bspline<point4d_t>  bspline_ops::reparametrize(const bspline<point4d_t>& spl,
                                             double t1, double t2);
//________________________________________________________
// method:split_into_bezier_patches_hard
template std::list<bspline<double>> bspline_ops::split_into_bezier_patches_hard(const bspline<double> &spl);
template std::list<bspline<point2d_t>> bspline_ops::split_into_bezier_patches_hard(const bspline<point2d_t> &spl);
template std::list<bspline<point3d_t>> bspline_ops::split_into_bezier_patches_hard(const bspline<point3d_t> &spl);
template std::list<bspline<point4d_t>> bspline_ops::split_into_bezier_patches_hard(const bspline<point4d_t> &spl);
//________________________________________________________
// method:split_into_bezier_patches
template std::list<bspline<double>> bspline_ops::split_into_bezier_patches(const bspline<double> &spl);
template std::list<bspline<point2d_t>> bspline_ops::split_into_bezier_patches(const bspline<point2d_t> &spl);
template std::list<bspline<point3d_t>> bspline_ops::split_into_bezier_patches(const bspline<point3d_t> &spl);
template std::list<bspline<point4d_t>> bspline_ops::split_into_bezier_patches(const bspline<point4d_t> &spl);
//________________________________________________________
// method:extract_curve
template bspline<double> bspline_ops::extract_curve(const bspline<double> &spl, double a, double b);
template bspline<point2d_t> bspline_ops::extract_curve(const bspline<point2d_t> &spl, double a, double b);
template bspline<point3d_t> bspline_ops::extract_curve(const bspline<point3d_t> &spl, double a, double b);
template bspline<point4d_t> bspline_ops::extract_curve(const bspline<point4d_t> &spl, double a, double b);
//________________________________________________________
// method:is_bezier
template bool  bspline_ops::is_bezier(const bspline<double>& spl);
template bool  bspline_ops::is_bezier(const bspline<point2d_t>& spl);
template bool  bspline_ops::is_bezier(const bspline<point3d_t>& spl);
template bool  bspline_ops::is_bezier(const bspline<point4d_t>& spl);
//________________________________________________________
// method:reverse_curve
template bspline<double>  bspline_ops::reverse_curve(bspline<double>& spl);
template bspline<point2d_t>  bspline_ops::reverse_curve(bspline<point2d_t>& spl);
template bspline<point3d_t>  bspline_ops::reverse_curve(bspline<point3d_t>& spl);
template bspline<point4d_t>  bspline_ops::reverse_curve(bspline<point4d_t>& spl);
//________________________________________________________
// method:foot_param
template double bspline_ops::foot_param(const bspline<double> &spl,
                        const  bspline<double>::point_t& p);
template double bspline_ops::foot_param(const bspline<point2d_t> &spl,
                        const  bspline<point2d_t>::point_t& p);
template double bspline_ops::foot_param(const bspline<point3d_t> &spl,
                        const  bspline<point3d_t>::point_t& p);
template double bspline_ops::foot_param(const bspline<point4d_t> &spl,
                        const  bspline<point4d_t>::point_t& p);
//________________________________________________________
// method:compute_box
template box<spline_traits<bspline<double>>::dim> bspline_ops::compute_box(const bspline<double> &spl);
template box<spline_traits<bspline<point2d_t>>::dim> bspline_ops::compute_box(const bspline<point2d_t> &spl);
template box<spline_traits<bspline<point3d_t>>::dim> bspline_ops::compute_box(const bspline<point3d_t> &spl);
template box<spline_traits<bspline<point4d_t>>::dim> bspline_ops::compute_box(const bspline<point4d_t> &spl);

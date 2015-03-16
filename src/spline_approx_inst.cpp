//________________________________________________________
// method:foot_param
template double bspline_ops::foot_param(const bspline<double> &spl,
                        const  bspline<double>::point_t& p) ;
template double bspline_ops::foot_param(const bspline<point2d_t> &spl,
                        const  bspline<point2d_t>::point_t& p) ;
template double bspline_ops::foot_param(const bspline<point3d_t> &spl,
                        const  bspline<point3d_t>::point_t& p) ;
template double bspline_ops::foot_param(const bspline<point4d_t> &spl,
                        const  bspline<point4d_t>::point_t& p) ;
template double bspline_ops::foot_param(const periodic_bspline<double> &spl,
                        const  periodic_bspline<double>::point_t& p) ;
template double bspline_ops::foot_param(const periodic_bspline<point2d_t> &spl,
                        const  periodic_bspline<point2d_t>::point_t& p) ;
template double bspline_ops::foot_param(const periodic_bspline<point3d_t> &spl,
                        const  periodic_bspline<point3d_t>::point_t& p) ;
template double bspline_ops::foot_param(const periodic_bspline<point4d_t> &spl,
                        const  periodic_bspline<point4d_t>::point_t& p) ;

//-*-mode:c++-*-
//Generated on: Tue Nov 17 14:00:19 2015. Do not edit
//________________________________________________________
// method:make_rbspline_from_circle
template rational_bspline< point2d_t, regular_tag > geom::make_rbspline_from_circle(const circle<point2d_t>& circ) ;
template rational_bspline< point3d_t, regular_tag > geom::make_rbspline_from_circle(const circle<point3d_t>& circ) ;
//________________________________________________________
// method:make_circle
template circle<point2d_t> geom::make_circle(const point2d_t& p1,
                  const point2d_t& p2,
                  const point2d_t& p3) ;
template circle<point3d_t> geom::make_circle(const point3d_t& p1,
                  const point3d_t& p2,
                  const point3d_t& p3) ;
//________________________________________________________
// method:foot_param
template double geom::foot_param(const circle<point2d_t> &c,
                 const point2d_t& p) ;
template double geom::foot_param(const circle<point3d_t> &c,
                 const point3d_t& p) ;

//-*-mode:c++-*-
//Generated on: Fri Aug 14 17:41:38 2015. Do not edit
//________________________________________________________
// method:insert_knots
template bspline<double> ops::insert_knots(const bspline<double>& crv,
                  const std::vector<double>& new_knots) ;
template bspline<point2d_t> ops::insert_knots(const bspline<point2d_t>& crv,
                  const std::vector<double>& new_knots) ;
template bspline<point3d_t> ops::insert_knots(const bspline<point3d_t>& crv,
                  const std::vector<double>& new_knots) ;
template bspline<point4d_t> ops::insert_knots(const bspline<point4d_t>& crv,
                  const std::vector<double>& new_knots) ;
template periodic_bspline<double> ops::insert_knots(const periodic_bspline<double>& crv,
                  const std::vector<double>& new_knots) ;
template periodic_bspline<point2d_t> ops::insert_knots(const periodic_bspline<point2d_t>& crv,
                  const std::vector<double>& new_knots) ;
template periodic_bspline<point3d_t> ops::insert_knots(const periodic_bspline<point3d_t>& crv,
                  const std::vector<double>& new_knots) ;
template periodic_bspline<point4d_t> ops::insert_knots(const periodic_bspline<point4d_t>& crv,
                  const std::vector<double>& new_knots) ;
template rational_bspline < point2d_t,regular_tag> ops::insert_knots(const rational_bspline < point2d_t,regular_tag>& crv,
                  const std::vector<double>& new_knots) ;
template rational_bspline < point3d_t,regular_tag> ops::insert_knots(const rational_bspline < point3d_t,regular_tag>& crv,
                  const std::vector<double>& new_knots) ;
template rational_bspline < double, regular_tag> ops::insert_knots(const rational_bspline < double, regular_tag>& crv,
                  const std::vector<double>& new_knots) ;
template rational_bspline < point2d_t,periodic_tag> ops::insert_knots(const rational_bspline < point2d_t,periodic_tag>& crv,
                  const std::vector<double>& new_knots) ;
template rational_bspline < point3d_t,periodic_tag> ops::insert_knots(const rational_bspline < point3d_t,periodic_tag>& crv,
                  const std::vector<double>& new_knots) ;
template rational_bspline < double,periodic_tag> ops::insert_knots(const rational_bspline < double,periodic_tag>& crv,
                  const std::vector<double>& new_knots) ;
//________________________________________________________
// method:insert_knot
template bspline<double> ops::insert_knot(const bspline<double>& crv,
                 double u) ;
template bspline<point2d_t> ops::insert_knot(const bspline<point2d_t>& crv,
                 double u) ;
template bspline<point3d_t> ops::insert_knot(const bspline<point3d_t>& crv,
                 double u) ;
template bspline<point4d_t> ops::insert_knot(const bspline<point4d_t>& crv,
                 double u) ;
template periodic_bspline<double> ops::insert_knot(const periodic_bspline<double>& crv,
                 double u) ;
template periodic_bspline<point2d_t> ops::insert_knot(const periodic_bspline<point2d_t>& crv,
                 double u) ;
template periodic_bspline<point3d_t> ops::insert_knot(const periodic_bspline<point3d_t>& crv,
                 double u) ;
template periodic_bspline<point4d_t> ops::insert_knot(const periodic_bspline<point4d_t>& crv,
                 double u) ;
template rational_bspline < point2d_t,regular_tag> ops::insert_knot(const rational_bspline < point2d_t,regular_tag>& crv,
                 double u) ;
template rational_bspline < point3d_t,regular_tag> ops::insert_knot(const rational_bspline < point3d_t,regular_tag>& crv,
                 double u) ;
template rational_bspline < double, regular_tag> ops::insert_knot(const rational_bspline < double, regular_tag>& crv,
                 double u) ;
template rational_bspline < point2d_t,periodic_tag> ops::insert_knot(const rational_bspline < point2d_t,periodic_tag>& crv,
                 double u) ;
template rational_bspline < point3d_t,periodic_tag> ops::insert_knot(const rational_bspline < point3d_t,periodic_tag>& crv,
                 double u) ;
template rational_bspline < double,periodic_tag> ops::insert_knot(const rational_bspline < double,periodic_tag>& crv,
                 double u) ;
//________________________________________________________
// method:match_knots
template bool ops::match_knots(bspline<double>& s1, bspline<double>& s2) ;
template bool ops::match_knots(bspline<point2d_t>& s1, bspline<point2d_t>& s2) ;
template bool ops::match_knots(bspline<point3d_t>& s1, bspline<point3d_t>& s2) ;
template bool ops::match_knots(bspline<point4d_t>& s1, bspline<point4d_t>& s2) ;
template bool ops::match_knots(periodic_bspline<double>& s1, periodic_bspline<double>& s2) ;
template bool ops::match_knots(periodic_bspline<point2d_t>& s1, periodic_bspline<point2d_t>& s2) ;
template bool ops::match_knots(periodic_bspline<point3d_t>& s1, periodic_bspline<point3d_t>& s2) ;
template bool ops::match_knots(periodic_bspline<point4d_t>& s1, periodic_bspline<point4d_t>& s2) ;
template bool ops::match_knots(rational_bspline < point2d_t,regular_tag>& s1, rational_bspline < point2d_t,regular_tag>& s2) ;
template bool ops::match_knots(rational_bspline < point3d_t,regular_tag>& s1, rational_bspline < point3d_t,regular_tag>& s2) ;
template bool ops::match_knots(rational_bspline < double, regular_tag>& s1, rational_bspline < double, regular_tag>& s2) ;
template bool ops::match_knots(rational_bspline < point2d_t,periodic_tag>& s1, rational_bspline < point2d_t,periodic_tag>& s2) ;
template bool ops::match_knots(rational_bspline < point3d_t,periodic_tag>& s1, rational_bspline < point3d_t,periodic_tag>& s2) ;
template bool ops::match_knots(rational_bspline < double,periodic_tag>& s1, rational_bspline < double,periodic_tag>& s2) ;

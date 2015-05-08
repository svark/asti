//-*-mode:c++-*-
//Generated on: Thu Mar 19 13:18:05 2015. Do not edit
//________________________________________________________
// class:bspline
template struct bspline<double>;
template struct bspline<point2d_t>;
template struct bspline<point3d_t>;
template struct bspline<point4d_t>;
//________________________________________________________
// method:blossom_eval,Point,KnotIter
template double  bspline<double>::blossom_eval(const double * f) const;
template double  bspline<double>::blossom_eval(std::vector<double>::const_iterator f) const;
template double  bspline<double>::blossom_eval(util::skip_ith_iter<std::vector<double>::const_iterator> f) const;
template point2d_t  bspline<point2d_t>::blossom_eval(const double * f) const;
template point2d_t  bspline<point2d_t>::blossom_eval(std::vector<double>::const_iterator f) const;
template point2d_t  bspline<point2d_t>::blossom_eval(util::skip_ith_iter<std::vector<double>::const_iterator> f) const;
template point3d_t  bspline<point3d_t>::blossom_eval(const double * f) const;
template point3d_t  bspline<point3d_t>::blossom_eval(std::vector<double>::const_iterator f) const;
template point3d_t  bspline<point3d_t>::blossom_eval(util::skip_ith_iter<std::vector<double>::const_iterator> f) const;
template point4d_t  bspline<point4d_t>::blossom_eval(const double * f) const;
template point4d_t  bspline<point4d_t>::blossom_eval(std::vector<double>::const_iterator f) const;
template point4d_t  bspline<point4d_t>::blossom_eval(util::skip_ith_iter<std::vector<double>::const_iterator> f) const;

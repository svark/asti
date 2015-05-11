//-*-mode:c++-*-
//Generated on: Fri May  8 16:20:08 2015. Do not edit
//________________________________________________________
// class:rational_bspline
template struct rational_bspline<point2d_t,regular_tag>;
template struct rational_bspline<point3d_t,regular_tag>;
template struct rational_bspline<double,regular_tag>;
template struct rational_bspline<point2d_t,periodic_tag>;
template struct rational_bspline<point3d_t,periodic_tag>;
template struct rational_bspline<double,periodic_tag>;
//________________________________________________________
// method:blossom_eval,Point,PTag,KnotIter
template  rational_bspline<point2d_t,regular_tag>::point_t rational_bspline<point2d_t,regular_tag>::blossom_eval(const double * us) const;
template  rational_bspline<point2d_t,regular_tag>::point_t rational_bspline<point2d_t,regular_tag>::blossom_eval(std::vector<double>::const_iterator us) const;
template  rational_bspline<point2d_t,regular_tag>::point_t rational_bspline<point2d_t,regular_tag>::blossom_eval(util::skip_ith_iter<std::vector<double>::const_iterator> us) const;
template  rational_bspline<point3d_t,regular_tag>::point_t rational_bspline<point3d_t,regular_tag>::blossom_eval(const double * us) const;
template  rational_bspline<point3d_t,regular_tag>::point_t rational_bspline<point3d_t,regular_tag>::blossom_eval(std::vector<double>::const_iterator us) const;
template  rational_bspline<point3d_t,regular_tag>::point_t rational_bspline<point3d_t,regular_tag>::blossom_eval(util::skip_ith_iter<std::vector<double>::const_iterator> us) const;
template  rational_bspline<double,regular_tag>::point_t rational_bspline<double,regular_tag>::blossom_eval(const double * us) const;
template  rational_bspline<double,regular_tag>::point_t rational_bspline<double,regular_tag>::blossom_eval(std::vector<double>::const_iterator us) const;
template  rational_bspline<double,regular_tag>::point_t rational_bspline<double,regular_tag>::blossom_eval(util::skip_ith_iter<std::vector<double>::const_iterator> us) const;
template  rational_bspline<point2d_t,periodic_tag>::point_t rational_bspline<point2d_t,periodic_tag>::blossom_eval(const double * us) const;
template  rational_bspline<point2d_t,periodic_tag>::point_t rational_bspline<point2d_t,periodic_tag>::blossom_eval(std::vector<double>::const_iterator us) const;
template  rational_bspline<point2d_t,periodic_tag>::point_t rational_bspline<point2d_t,periodic_tag>::blossom_eval(util::skip_ith_iter<std::vector<double>::const_iterator> us) const;
template  rational_bspline<point3d_t,periodic_tag>::point_t rational_bspline<point3d_t,periodic_tag>::blossom_eval(const double * us) const;
template  rational_bspline<point3d_t,periodic_tag>::point_t rational_bspline<point3d_t,periodic_tag>::blossom_eval(std::vector<double>::const_iterator us) const;
template  rational_bspline<point3d_t,periodic_tag>::point_t rational_bspline<point3d_t,periodic_tag>::blossom_eval(util::skip_ith_iter<std::vector<double>::const_iterator> us) const;
template  rational_bspline<double,periodic_tag>::point_t rational_bspline<double,periodic_tag>::blossom_eval(const double * us) const;
template  rational_bspline<double,periodic_tag>::point_t rational_bspline<double,periodic_tag>::blossom_eval(std::vector<double>::const_iterator us) const;
template  rational_bspline<double,periodic_tag>::point_t rational_bspline<double,periodic_tag>::blossom_eval(util::skip_ith_iter<std::vector<double>::const_iterator> us) const;

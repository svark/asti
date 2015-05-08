//-*-mode:c++-*-
//Generated on: Fri Mar 27 19:34:16 2015. Do not edit
//________________________________________________________
// class:rational_bspline
template struct rational_bspline<bspline<point2d_t>>;
template struct rational_bspline<bspline<point3d_t>>;
template struct rational_bspline<bspline<point4d_t>>;
template struct rational_bspline<periodic_bspline<point2d_t>>;
template struct rational_bspline<periodic_bspline<point3d_t>>;
template struct rational_bspline<periodic_bspline<point4d_t>>;
//________________________________________________________
// method:blossom_eval,SplineType,KnotIter
template  rational_bspline<bspline<point2d_t>>::point_t rational_bspline<bspline<point2d_t>>::blossom_eval(const double * us) const;
template  rational_bspline<bspline<point2d_t>>::point_t rational_bspline<bspline<point2d_t>>::blossom_eval(std::vector<double>::const_iterator us) const;
template  rational_bspline<bspline<point2d_t>>::point_t rational_bspline<bspline<point2d_t>>::blossom_eval(util::skip_ith_iter<std::vector<double>::const_iterator> us) const;
template  rational_bspline<bspline<point3d_t>>::point_t rational_bspline<bspline<point3d_t>>::blossom_eval(const double * us) const;
template  rational_bspline<bspline<point3d_t>>::point_t rational_bspline<bspline<point3d_t>>::blossom_eval(std::vector<double>::const_iterator us) const;
template  rational_bspline<bspline<point3d_t>>::point_t rational_bspline<bspline<point3d_t>>::blossom_eval(util::skip_ith_iter<std::vector<double>::const_iterator> us) const;
template  rational_bspline<bspline<point4d_t>>::point_t rational_bspline<bspline<point4d_t>>::blossom_eval(const double * us) const;
template  rational_bspline<bspline<point4d_t>>::point_t rational_bspline<bspline<point4d_t>>::blossom_eval(std::vector<double>::const_iterator us) const;
template  rational_bspline<bspline<point4d_t>>::point_t rational_bspline<bspline<point4d_t>>::blossom_eval(util::skip_ith_iter<std::vector<double>::const_iterator> us) const;
template  rational_bspline<periodic_bspline<point2d_t>>::point_t rational_bspline<periodic_bspline<point2d_t>>::blossom_eval(const double * us) const;
template  rational_bspline<periodic_bspline<point2d_t>>::point_t rational_bspline<periodic_bspline<point2d_t>>::blossom_eval(std::vector<double>::const_iterator us) const;
template  rational_bspline<periodic_bspline<point2d_t>>::point_t rational_bspline<periodic_bspline<point2d_t>>::blossom_eval(util::skip_ith_iter<std::vector<double>::const_iterator> us) const;
template  rational_bspline<periodic_bspline<point3d_t>>::point_t rational_bspline<periodic_bspline<point3d_t>>::blossom_eval(const double * us) const;
template  rational_bspline<periodic_bspline<point3d_t>>::point_t rational_bspline<periodic_bspline<point3d_t>>::blossom_eval(std::vector<double>::const_iterator us) const;
template  rational_bspline<periodic_bspline<point3d_t>>::point_t rational_bspline<periodic_bspline<point3d_t>>::blossom_eval(util::skip_ith_iter<std::vector<double>::const_iterator> us) const;
template  rational_bspline<periodic_bspline<point4d_t>>::point_t rational_bspline<periodic_bspline<point4d_t>>::blossom_eval(const double * us) const;
template  rational_bspline<periodic_bspline<point4d_t>>::point_t rational_bspline<periodic_bspline<point4d_t>>::blossom_eval(std::vector<double>::const_iterator us) const;
template  rational_bspline<periodic_bspline<point4d_t>>::point_t rational_bspline<periodic_bspline<point4d_t>>::blossom_eval(util::skip_ith_iter<std::vector<double>::const_iterator> us) const;

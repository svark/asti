//-*-mode:c++-*-
//Generated on: Tue Mar 10 18:02:19 2015. Do not edit
//________________________________________________________
// class:rational_bspline
template class rational_bspline<bspline<point2d_t>>;
template class rational_bspline<bspline<point3d_t>>;
template class rational_bspline<bspline<point4d_t>>;
template class rational_bspline<periodic_bspline<point2d_t>>;
template class rational_bspline<periodic_bspline<point3d_t>>;
template class rational_bspline<periodic_bspline<point4d_t>>;
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

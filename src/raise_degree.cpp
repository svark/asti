#include "skip_idx_iter.hpp"
#include "raise_degree.hpp"
#include "rational_bspline_cons.hpp"
namespace geom {

//{{{ --(@* "raise the degree by 1")
// see:(@file :file-name "media/raise_degree.png" :to "./media/raise_degree.png" :display "raise_degree")
template <class SplineType>
SplineType ops::raise_degree(const SplineType& spl)
{
    typedef typename SplineType::knots_t knots_t;
    //typedef typename SplineType::point_t Point;
    typedef typename SplineType::vector_t vector_t;

    knots_t  new_knots;
    new_knots.reserve(spl.knots().size());
    auto b = spl.knots().cbegin();
    auto e = spl.knots().cend();

    for(auto f = b; f != e; ++f) {
        new_knots.push_back(*f);
        new_knots.push_back(*f);
    }

    size_t num_new_knots = new_knots.size();
    int p = spl.degree() + 1;
    size_t num_new_cpts = num_new_knots - p - 1;

    using util::skip_ith_iter;
    SplineType::cpts_t new_cpts;
    new_cpts.reserve(num_new_cpts);
    for(size_t i = 0; i < num_new_cpts; ++i)  {
        vector_t cv(0.0);
        for(int j = 0; j < p; ++j) {
            skip_ith_iter<decltype(new_knots.cbegin())>
                iter( j, new_knots.cbegin() + (i + 1));
            cv += make_vec(spl.blossom_eval(iter)) ;
        }
        cv *= (1.0/p);
        new_cpts.push_back(make_pt(cv));
    }
    return SplineType(std::move(new_cpts),
                      std::move(new_knots), p).translate(spl.base_point());
}

template <class SplineCurve>
rational_bspline <SplineCurve>
raise_degree(const rational_bspline<SplineCurve>
           & crv, double u)
{
    return make_rbspline(raise_degree(crv.spline()));
}
//}}}
}
//{{{  instantiation scripts

/*
  Local Variables:
  eval:(load-file "./scripts/temp.el")
  eval:(setq methods (list "raise_degree"
  ))
  eval:(setq spltypes (list "bspline<double>"
  "bspline<point2d_t>"
  "bspline<point3d_t>"
  "bspline<point4d_t>"
  "periodic_bspline<double>"
  "periodic_bspline<point2d_t>"
  "periodic_bspline<point3d_t>"
  "periodic_bspline<point4d_t>"
  "rational_bspline < bspline<point2d_t>>"
  "rational_bspline < bspline<point3d_t>>"
  "rational_bspline < bspline<point4d_t>>"
  "rational_bspline < periodic_bspline<point2d_t>>"
  "rational_bspline < periodic_bspline<point3d_t>>"
  "rational_bspline < periodic_bspline<point4d_t>>"
  ))
  eval:(instantiate-templates "raise_degree" "ops" (list ) methods spltypes )
  End:
// dump all explicitly instantiated templates below
*/
//}}}

//{{{  instantiation
#include "bspline.hpp"
#include "periodic_bspline.hpp"
#include "point.hpp"
namespace geom {
#include "raise_degree_inst.inl"
}
//}}}

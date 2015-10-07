//-*- mode:c++ -*-
#include "skip_idx_iter.hpp"
#include "raise_degree.hpp"
#include "bspline_x_cons.hpp"
#include <algorithm>
#include "smat.hpp"

namespace geom {

//{{{ --(@* "raise the degree by 1")
// media/raise_degree.png
// media/raise_degree_1.png
namespace impl
{
template <class SplineType>
SplineType raise_degree_helper(const SplineType& spl, polynomial_tag, regular_tag)
{
    typedef typename SplineType::knots_t knots_t;
    //typedef typename SplineType::point_t Point;
    typedef typename SplineType::vector_t vector_t;

    knots_t  new_knots;
    new_knots.reserve(2*spl.knots().size());


    auto comp = [](double u, double v)
        {
            return  u < v - tol::param_tol/2;
        };
    auto s = clamp_end(clamp_start(spl));
    auto b = s.knots().cbegin();
    auto e = s.knots().cend();

    for(auto f = b; f != e;) {
        auto l = std::upper_bound(f,e,*f, comp);
        std::copy(f,l, std::back_inserter(new_knots));
        new_knots.push_back(*f);
        f = l;
    }

    size_t num_new_knots = new_knots.size();
    int p = spl.degree() + 1;
    size_t num_new_cpts = num_new_knots - p - 1;

    using util::skip_ith_iter;
    typename SplineType::cpts_t new_cpts;
    new_cpts.reserve(num_new_cpts);
    for(size_t i = 0; i < num_new_cpts; ++i)  {
        vector_t cv(0.0);
        for(int j = 0; j < p; ++j) {
            skip_ith_iter<decltype(new_knots.cbegin())>
                iter( j, new_knots.cbegin() + (i + 1));
            std::vector<double> tmp(p-1);
            std::copy_n(iter, p-1,tmp.begin() );
            cv += make_vec(spl.blossom_eval(tmp.cbegin())) ;
        }
        cv *= (1.0/p);
        new_cpts.push_back(make_pt(cv));
    }
    typedef spline_traits<SplineType> spl_traits;
    return make_bspline(std::move(new_cpts),
                        std::move(new_knots), p,
                        typename spl_traits::ptag(),
                        typename spl_traits::rtag()
        );

}

template <class SplineType>
SplineType raise_degree_helper(const SplineType& spl, polynomial_tag, periodic_tag)
{
    int p = spl.degree();
    std::vector<double> sknots,eknots;
    sknots.reserve(p + 2); eknots.reserve(p + 2);
    auto const &tsb =  spl.knots().begin();
    sknots.assign(tsb, tsb + p + 1);
    sknots.push_back(tsb[p]);

    auto const &tse =  spl.knots().end() - p - 1 ;
    eknots.push_back(tse[0]);
    eknots.insert(eknots.end(),tse, tse + p + 1);

    auto bs = clamp_start(spl);
    clamp_end(bs).swap(bs);
    raise_degree_helper(bs, polynomial_tag(), regular_tag()).swap(bs);
    rebase_at_start(bs,sknots.begin()).swap(bs);

    rebase_at_end(bs, eknots.begin()).swap(bs);
    return bs;
}

template <class SplineType>
SplineType raise_degree_helper(const SplineType& spl, rational_tag, regular_tag)
{
    return make_rbspline( raise_degree_helper(spl.spline(), polynomial_tag(), regular_tag()) );
}

template <class SplineType>
SplineType raise_degree_helper(const SplineType& spl, rational_tag, periodic_tag)
{
    return make_rbspline( raise_degree_helper(spl.spline(), polynomial_tag(), periodic_tag()) );
}
}
template <class SplineType>
SplineType ops::raise_degree(const SplineType& spl)
{
    return impl::raise_degree_helper(spl, typename spline_traits<SplineType>::rtag(),
                                     typename spline_traits<SplineType>::ptag()
        );
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
  "rational_bspline < point2d_t,regular_tag>"
  "rational_bspline < point3d_t,regular_tag>"
  "rational_bspline < double, regular_tag>"
  "rational_bspline < point2d_t,periodic_tag>"
  "rational_bspline < point3d_t,periodic_tag>"
  "rational_bspline < double,periodic_tag>"
  ))
  eval:(instantiate-templates "raise_degree" "ops" (list ) (product methods
  spltypes ))
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

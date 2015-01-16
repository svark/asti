//--*-mode:c++-mode-*-
//{{{ --includes
#include "stdafx.h"
#include <iostream>
#include <map>
#include <cmath>

#include "bspline.hpp"
#include "rmat.hpp"
#include "tol.hpp"
#include "util.hpp"
#include <utility>
#include <memory>
//}}}
namespace geom {

  //{{{ --(@* "evaluate curve at param @u")
  template <class Point>
  Point bspline<Point>::eval(double u) const
  {
    if(tol::eq(u,t.back(),tol::param_tol))
      return cpts.back();

    rmat<Point> m(cpts,t, deg);
    return m.eval(u) + origin;
  }
  //}}}

  //{{{ -- (@* "evaluate derivative of the curve at the param @u")
  // throws a spline_exception if u is not in parameter range of this curve
  template <class Point>
  typename bspline<Point>::vector_t
  bspline<Point>::eval_derivative(int numDer, double u) const
  {
    rmat<Point> m(cpts, t, deg);
    assert(numDer >= 0);
    return  make_vec( m.eval_derivative(numDer, u) ) ;
  }

  //}}}

  //{{{ --(@* "parameter range")
  template <class Point>
  std::pair<double,double>
  bspline<Point>::param_range() const
  {
    return std::make_pair(t.front(), t.back() );
  }
  //}}}

  //{{{ --(@* "constructors")

  template <class Point>
  bspline<Point>::bspline(cpts_t&& pts,
                          knots_t &&ks, int degree_):
    cpts(std::move(pts)),t(std::move(ks)),deg(degree_),
    origin(vector_t(0.0))
  {
  }

  template <class Point>
  bspline<Point>::bspline(const cpts_t& pts,
                          const knots_t &ks, int degree_):
    cpts(pts),t(ks),deg(degree_),
    origin(vector_t(0.0))
  {
  }

  template <class Point>
  bspline<Point>::bspline(const bspline<Point> &other):
    cpts(other.cpts),t(other.t),deg(other.deg),
    origin(other.origin)
  {
  }

  template <class Point>
  bspline<Point>::bspline(bspline<Point>&& other):
    cpts(std::move(other.cpts)),t(std::move(other.t)),deg(other.deg),
    origin(other.origin)
  {

  }

  template <class Point>
  void bspline<Point>::optimize()
  {
    auto center = make_vec(centroid( cpts.begin(), cpts.end() ) );
    std::transform(cpts.begin(), cpts.end(), cpts.begin(),
                   [&center](Point& p){ return p -= center;} );
    origin += center;
  }
  //}}}

  //{{{ --(@*"evaluate bspline blossom f[0],f[1],...,f[deg-1]")
  //see
  //(@file :file-name "./blossom1.png" :to "./blossom1.png" :display "blossom1")
  //(@file :file-name "./blossom2.png" :to "./blossom2.png" :display "blossom2")
  template <class Point>
  template <class knot_iter>
  Point
  bspline<Point>::blossom_eval(knot_iter f)
  {
    int p = degree();
    std::unique_ptr<point_t[]> cache(new point_t[p+1]);
    rmat<Point> m(cpts,t,deg);
    size_t nu = m.locate_nu(f[0]);
    for(int j = 0;j < p + 1; ++j) {
      cache[j] = control_points()[nu - p + j];
    }
    m.spline_compute(f, nu, cache.get());
    return cache[0] + origin;
  }
  //}}}


}

using geom::bspline;

#include "point.hpp"

using geom::pt_t;


template struct bspline<double>;
template struct bspline<pt_t<2> >;
template struct bspline<pt_t<3> >;
template struct bspline<pt_t<4> >;

template pt_t<2> bspline<pt_t<2>>::blossom_eval( const double *);
template pt_t<3> bspline<pt_t<3>>::blossom_eval( const double *);
template pt_t<4> bspline<pt_t<4>>::blossom_eval( const double *);
template double bspline<double>::blossom_eval( const double *);

template pt_t<2> bspline<pt_t<2>>::blossom_eval( std::vector<double>::const_iterator);
template pt_t<3> bspline<pt_t<3>>::blossom_eval( std::vector<double>::const_iterator);
template pt_t<4> bspline<pt_t<4>>::blossom_eval( std::vector<double>::const_iterator);
template double bspline<double>::blossom_eval( std::vector<double>::const_iterator);

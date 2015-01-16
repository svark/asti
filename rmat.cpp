//-*- mode:c++ -*-
#include <stdafx.h>
#include <cassert>
#include "geom_exception.hpp"
//#include <algorithm>
#include "tol.hpp"
#include "rmat.hpp"
#include <algorithm>
#include <memory>
#include <math.h>
namespace geom {

  rmat_base::knots_t::const_iterator
  rmat_base::locate(double u) const
  {
    auto it = std::upper_bound( t.cbegin(),
                                t.cend(), u);
    if (it != t.cend() && it != t.cbegin() )
      return std::prev(it);

    if( it == t.cend() && fabs(u - t.back() ) < tol::param_tol ) {
      it = std::lower_bound(t.cbegin() + deg, t.cend(),
                                        t.back());
	}
	else if( it == t.cbegin() && fabs(u - t.front()) < tol::param_tol)
	{
		it = std::upper_bound( t.cbegin(),
                               t.cend(), t.front());
	}
	if( it != t.cend() && it != t.cbegin() )
		return std::prev(it);

    throw geom_exception(knot_not_in_range_error);
  }

  double
  rmat_base::der_n(size_t idx,
                   int numDer,
                   double u) const
  {
    int size = degree();
    auto nu = locate_nu(u);
    if(idx < nu - size || idx > nu )
      return 0;
    std::unique_ptr<double[]> cache(new double[size+1]);
    std::fill(cache.get(), cache.get() + size +  1,  0.0);

    cache[idx-nu+size] = 1;
    spline_compute(nu, u, numDer, cache.get());
    return cache[0];
  }

  template <class T>
  void
  rmat_base::spline_compute(size_t nu,
                            double u,
                            int numDer,
                            T  cache) const
  {
    int size = deg;
    auto t_ =  t.cbegin() + nu;
    size_t fac = 1;

	if(numDer > deg ) 
		cache[0] = std::remove_reference<decltype(cache[0])>::type(0.0);
    for(int sz = size; sz > size - numDer; --sz) {
      fac *= sz;
      for(size_t j = 1; j <= sz; ++j) {
        double r = t_[j] - t_[j - sz];
        if(r < tol::param_tol)
          throw geom_exception(knot_not_in_range_error_der);

        double lambda =  1 / r;
        cache[j - 1] = dlerp(lambda, cache[j - 1], cache[j]);
      }
    }

    for(int sz = size - numDer; sz > 0; --sz) {
      for(size_t j = 1; j <= sz; ++j) {
        double r = t_[j] - t_[j - sz];
        if(r < tol::param_tol)
          continue;

        assert(u >= t_[j - sz]);
        double lambda =  (u  - t_[j - sz]) / r;
        cache[j - 1] = lerp(lambda, cache[j - 1], cache[j] );
      }
    }
    scale(cache[0],double(fac));
  }



  template <class Point>
  Point
  rmat<Point>::eval_derivative(int numDer, double u) const
  {
    size_t nu = locate_nu(u);

    int size = degree();

    std::unique_ptr<point_t[]> cache(new point_t[size+1]);

    for(size_t j = 0; j < size + 1; ++j)
      cache[j] = wcpt(j + nu - size);

    spline_compute(nu, u, numDer, cache.get());

    return cache[0];
  }

  template <class Point>
  template <class ParamIter,class PointOutIter>
  void
  rmat<Point>::eval_derivatives(int numDer,
                                ParamIter us,
                                ParamIter end,
                                PointOutIter out)
  {
    int size = degree();

    std::unique_ptr<point_t[]> cache(new point_t[size+1]);

    size_t nu = -1;
    for( ;us != end; ++us,++out) {

      nu = locate_nu(u,nu);

      for(int j = 0; j < size + 1; ++j)
        cache[j] = wcpt(j + nu - size);

      spline_compute(nu, u, numDer, cache.get());

      *out =  cache[0];
    }

  }

  template <class Point>
  typename rmat<Point>::cpts_t
  rmat<Point>::insert_knots(const std::vector<double>& taus)
  {
    std::vector<point_t> newcpts;
    int size = degree();
    newcpts.reserve( taus.size() - size - 1 );

    for(size_t i = 0; i < taus.size() - size - 1; ++i) {

      size_t nu = locate_nu(taus[i]);

      std::unique_ptr<point_t[]> cache(new point_t[size + 1]);

      for(size_t j = 0; j < size + 1; ++j)
        cache[j] = wcpt(j + nu - size);

      spline_compute(taus.cbegin() + i, nu, cache.get());

      newcpts.push_back(cache[0]);
    }
    return newcpts;
  }
}
#include "point.hpp"

template  struct geom::rmat<double>;
template  struct geom::rmat<geom::pt_t<2>>;
template  struct geom::rmat<geom::pt_t<3>>;
template  struct geom::rmat<geom::pt_t<4>>;

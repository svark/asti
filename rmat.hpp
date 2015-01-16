//-*- mode:c++ -*-
#ifndef RMAT_HPP
#define RMAT_HPP
#include <vector>

namespace geom {

  struct rmat_base {
    typedef std::vector<double> knots_t;

    template <class KnotIter>
    rmat_base (KnotIter b, KnotIter e, int deg_)
      :t(b,e), deg(deg_) {
    }
    rmat_base ( const knots_t& t_, int deg_)
      :t(t_), deg(deg_) {
    }

    // find knot range [t[mu],t[mu+1]) in which @u belongs
    knots_t::const_iterator
    locate(double u) const;

    // find degree of the basis polynomials
    int
    degree() const {
      return deg;
    }

    const knots_t &
    knots() const {
      return t;
    }

    double
    coeff(size_t idx,  double u) const {
      return der_n(idx, 0, u);
    }
    double
    der(size_t idx,  double u) const {
      return der_n(idx, 1, u);
    }
    template <class T>
    void spline_compute(size_t nu,
                        double u,
                        int numDer,
                        T cache) const;

    template <class KnotIter, class T>
    void
    spline_compute(KnotIter us,
                   size_t nu, //location of us[0] in knot seq t
                   T cache) const

    {
      int size = deg;
      auto t_ =  t.cbegin() + nu;
      assert(nu >= size);

      for(int sz = size; sz > 0; --sz) {
        for(int j = 1; j <= sz; ++j) {

          double r = t_[j] - t_[j - sz];

          if(r < tol::param_tol)
            continue;

          double lambda =  (us[sz]  - t_[j - sz]) / r;

          if(lambda < 0 || lambda > 1)
            continue;

          cache[j - 1] = lerp(lambda,
                              cache[j - 1], cache[j] );
        }
      }
    }

    double
    der_n(size_t idx,int numDer, double u) const;

    size_t
    locate_nu(double u) const {
      auto b = knots().cbegin();
      auto nu = std::distance(b,
                              locate(u));
      return nu;
    }
    size_t
    locate_nu(double u, size_t nu_guess) const {
      auto b = knots().cbegin();
      if( b[nu_guess] < u && u < b[nu_guess+1] )
        return nu_guess;
      return locate_nu(u);
    }
  private:
    const knots_t& t;
    int deg;
  };

  template<class Point>
  struct rmat : public rmat_base {

    typedef Point point_t;
    typedef std::vector<point_t> cpts_t;

    rmat(const cpts_t& cpts_,const knots_t& t,int deg)
      : rmat_base(t,deg),cpts(cpts_)
    {
    }

    // evaluate the bspline curve at parameter u
    point_t
    eval(double u) const {
      return eval_derivative(0,u);
    }

    point_t
    eval_derivative(int numDer, double u) const;

    template <class ParamIter,class PointIter>
    void
    eval_pts(ParamIter f, ParamIter l, PointIter out) const {
      return eval_derivatives(0,f,l, out);
    }

    template <class ParamIter,class PointOutIter>
    void eval_derivatives(int numDer,
                          ParamIter us, ParamIter end,
                          PointOutIter out);

    std::vector<point_t>
    insert_knots(const std::vector<double>& taus);

    const cpts_t& control_pts() const {
      return cpts;
    }
    const point_t& wcpt(size_t i) const {
      return cpts[i];
    }
    const cpts_t& cpts;
  };
}
#endif //RMAT_HPP

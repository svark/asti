//-*- mode:c++ -*-
#ifndef ASTI_RMAT_HPP
#define ASTI_RMAT_HPP
#include <vector>
#include <memory>
#include "rmat_explicit.hpp"
#include <Eigen/Dense>
#include "point_fwd.hpp"
#include "type_utils.hpp"
namespace geom {
template <class KnotIter>
struct rmat_base {

    rmat_base(KnotIter b_,
              KnotIter e_,
              int deg_)
        :t(b_), e(e_),
         deg(deg_)
    {
    }

    //  knot range [t[mu],t[mu+1]) in which @u belongs
    KnotIter
    locate(double u) const;

    double front() const {  return  t[0]; } ;
    double back() const { return   e[-1] ;};

    size_t start_mult() ;
    size_t end_mult();
    size_t mult(double u);

    void swap(rmat_base&o)
    {
        std::swap(t,o.t);
        std::swap(e,o.e);
        std::swap(deg,o.deg);
    }
    // degree of the basis polynomials
    int
    degree() const {
        return deg;
    }

    double
    coeff(size_t idx,  double u) const {
        return der_n(idx, 0, u);
    }

    double
    der(size_t idx,  double u) const {
        return der_n(idx, 1, u);
    }

    // bspline derivative coefficients optimized for space
    double der_n(size_t idx, int derOrder, double u) const;

    // get the array of deg + 1 bspline basis functions at u
    //.  ./media/basis_comp.png
    // returns the vector b = (B_{\mu - p, p}(u), \ldots, B_{\mu, p}(x)), mu =
    // int such that u\in [t_\mu, t_\mu + 1) and p is the degree
    //also  returned is the index in the knot array -> \mu.
    std::tuple<std::vector<double>, size_t>
    get_basis(double u,int derOrder = 0) const;

    // knot insertion matrix of size:(l - f) * (e - t)
    Eigen::MatrixXd
    insertion_matrix(KnotIter f, KnotIter l) const;

    // locate the greatest index nu such that t[nu] <= u, except when
    // u == back() in which case return greatest nu such that
    // t[nu] <= u - tol
    size_t locate_nu(double u) const {
        auto nu = std::distance(t, locate(u));
        return nu;
    }
    // same as above but with a guess for nu to aid computation
    size_t
    locate_nu(double u, size_t nu_guess) const;

    template <class PointIter>
    bool eval(size_t nu, double u,
              int derOrder, PointIter cache) const;

    //leaner implementation of get_basis
    std::pair<std::vector<double>,size_t> basis(double u);

    template <class KnotIterU, class PointIter>
    void blossom_eval(KnotIterU us,
                      size_t nu, //location of us[0] in
                      //knot seq t
                      PointIter cache) const;

	template <class PointIter>
	bool der_eval(size_t nu, int derOrder, PointIter cache) const;

protected:
    KnotIter t,e;
    int deg;
};


template <class KnotIter>
template <class PointIter>
bool
rmat_base<KnotIter>::der_eval(size_t nu,
                              int derOrder,
                              PointIter cache) const
{
    int size = deg;

    auto t_ =  t + nu;
    size_t fac = 1;

    typedef RAWTYPE(cache[0]) point_t;
    if(derOrder > size )
        cache[0] = point_t(0.0);
	else
    for(int sz = size; sz > size - derOrder; --sz) {
        fac *= sz;
        for(int j = 1; j <= sz; ++j) {
            double d = t_[j] - t_[j - sz];
            if(tol::small(d, tol::param_tol)) {
                set_quiet_NaN(cache[0]);
                return false;
            }

            double lambda =  1 / d;
            cache[j - 1] = dlerp(lambda, cache[j - 1], cache[j]);
        }
    }
    for(int i = 0 ; i <= size - derOrder; ++i)
    {
        cache[i] *= double(fac);
    }
    return true;
}

template <class KnotIter>
template <class PointIter>
bool
rmat_base<KnotIter>::eval(size_t nu,
                          double u,
                          int derOrder,
                          PointIter cache) const
{
    assert(t[nu] <= u || u < t[deg] );
    int size = deg;
    auto t_ =  t + nu;

    if(!der_eval(nu,derOrder,cache))
        return false;

    for(int sz = size - derOrder; sz > 0; --sz) {
        for(int j = 1; j <= sz; ++j) {
            double d = t_[j] - t_[j - sz];
            double lambda =  sdiv(u  - t_[j - sz], d);
            cache[j - 1] = lerp(lambda, cache[j - 1], cache[j] );
        }
    }
    return true;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <class KnotIter>
template <class KnotIterU, class PointIter>
void
rmat_base<KnotIter>::blossom_eval( 
    KnotIterU us,  // increasing array of d knots
    size_t nu,     // location of us[0] in knot seq t
    PointIter cache) const
{
    assert(t[nu] <= *us || *us < t[deg] );
    int size = deg;
    auto t_ =  t + nu;
    assert(nu >= size_t(size));
    for(int sz = size; sz > 0; --sz) {
        for(int j = 1; j <= sz; ++j) {
            double d = t_[j] - t_[j - sz];
            double lambda =  sdiv(us[sz-1]  - t_[j - sz], d);
            cache[j - 1] = lerp(lambda, cache[j - 1], cache[j] );
        }
    }
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
struct rmat_base_vd: rmat_base < std::vector<double>::const_iterator >
{
    typedef std::vector<double> knots_t;
    rmat_base_vd(const knots_t & t, int deg):
        rmat_base(t.cbegin(), t.cend(), deg){}

};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template<class Point>
struct rmat : public rmat_base_vd {
    typedef Point point_t;
    typedef decltype(mk_stdvec(point_t())) cpts_t;

    rmat(const cpts_t& cpts_,const knots_t& t,int deg)
        : rmat_base_vd(t, deg), cpts(cpts_)
    {
    }

    point_t
    eval(double u) const {
        return eval_derivative(0,u);
    }

    point_t
    eval_derivative(int derOrder, double u) const;

    template <class ParamIter,class PointIter>
    void
    eval(ParamIter f, ParamIter l, PointIter out) const {
        return eval_derivatives(0,f,l, out);
    }

    template <class ParamIter,class PointOutIter>
    void
    eval_derivatives(int derOrder,
                     ParamIter us, ParamIter end,
                     PointOutIter out) const
    {
        int size = degree();
        std::unique_ptr<point_t[]> cache(new point_t[size+1]);

        size_t nu = -1;
        for( ;us != end; ++us,++out) {

            nu = locate_nu( * us,nu);

            for(int j = 0; j < size + 1; ++j)
                cache[j] = control_pt(j + nu - size);

            eval(nu,* us, derOrder, cache.get());

            *out =  cache[0];
        }
    }

    template <class PointOutIter>
    void
    eval_derivatives(int derOrder,
                     double u,
                     PointOutIter out)
    {
        int size = degree();
        std::unique_ptr<point_t[]> cache(new point_t[size+1]);
        size_t nu = locate_nu(u);

        for(int i =  0; i <= derOrder;++i,++out) {

            for(int j = 0; j < size + 1; ++j)
                cache[j] = control_pt(j + nu - size);

            rmat_base_vd::eval(nu, u, i, cache.get());

            *out = cache[0];
        }
    }

    cpts_t
    insert_knots(const std::vector<double>& taus);

    const cpts_t& control_pts() const {
        return cpts;
    }
    point_t control_pt(size_t i) const {
        return cpts[i];
    }
    const cpts_t& cpts;
};
}
#endif //RMAT_HPP

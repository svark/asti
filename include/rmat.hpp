//-*- mode:c++ -*-
#ifndef ASTI_RMAT_HPP
#define ASTI_RMAT_HPP
#include <vector>
#include <memory>
#include "rmat_explicit.hpp"
#include <Eigen/Dense>
#include "point.hpp"
#include "geom_exception.hpp"
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
    // function is optimized for parallelism
    std::vector<double> der_coeffs_par(int derOrder,
                                       double u) const;

    // deg + 1, bspline coefficients at u, optimized for parallelism
    std::vector<double> coeffs_par(double u) const;

    template<class Point>
    struct accumulator  {
        accumulator(double u, int derOrder):base(0.0), k(0)
        {
            size_t nu = locate_nu(u);
            assert(derOrder >= 0);
            for(int j = 1;j < deg - derOrder; ++j) {
                basis_cache *= rmat_explicit<KnotIter>(t + nu, j, u);
            }
            for(int j = deg - derOrder;j <= deg; ++j) {
                basis_cache *= der_rmat_explicit<KnotIter>(t + nu, j, u);
            }
        }

        accumulator(accumulator&& other)
            :basis_cache(other.basis_cache)
        {
            k = 0;
            base = Point(0.0);
        }

        void prod(const Point& pt)
        {
            base = axpy(basis_cache.get(k++), pt, base);
        }

        Point get() const
        {
            return base;
        }

        void swap(accumulator &other)
        {
            basis_cache.swap(other.basis_cache);
            std::swap(k, other.k);
            std::swap(base, other.base);
        }
    private:
        mult_rmat basis_cache;
        int k;
        Point base;
    };

    template<class Point>
    accumulator < Point >
    get_accumulator(double u, int derOrder = 0)
    {
        return accumulator < Point > (u, derOrder);
    }

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

    template <class T>
    void spline_compute(size_t nu, double u, int derOrder, T cache) const;

    //.  ./media/basis_comp.png
    // returns the vector b = (B_{\mu - p, p}(u), \ldots, B_{\mu, p}(x)), mu =
    // int such that u\in [t_\mu, t_\mu + 1) and p is the degree
    std::vector<double> basis(double u);

    template <class KnotIterU, class PointIter>
    void spline_compute(KnotIterU us, size_t nu, //location of us[0] in
                                                //knot seq t
                        PointIter cache) const;

protected:
    KnotIter t,e;
    int deg;
};

template <class KnotIter>
template <class PointIter>
void
rmat_base<KnotIter>::spline_compute(size_t nu,
                        double u,
                        int derOrder,
                        PointIter cache) const
{
    int size = deg;
    auto t_ =  t + nu;
    size_t fac = 1;

    typedef RAWTYPE(cache[0]) point_t;
    if(derOrder > size )
        cache[0] = point_t(0.0);

    for(int sz = size; sz > size - derOrder; --sz) {
        fac *= sz;
        for(int j = 1; j <= sz; ++j) {
            double d = t_[j] - t_[j - sz];
            if(d < tol::param_tol)
                throw geom_exception(knot_not_in_range_error_der);

            double lambda =  1 / d;
            cache[j - 1] = dlerp(lambda, cache[j - 1], cache[j]);
        }
    }

    for(int sz = size - derOrder; sz > 0; --sz) {
        for(int j = 1; j <= sz; ++j) {
            double d = t_[j] - t_[j - sz];
            double lambda =  sdiv(u  - t_[j - sz], d);
            cache[j - 1] = lerp(lambda, cache[j - 1], cache[j] );
        }
    }
    cache[0] *= double(fac);
}

template <class KnotIter>
template <class KnotIterU, class PointIter>
void
rmat_base<KnotIter>::spline_compute(
    KnotIterU us,
    size_t nu, //location of us[0] in knot seq t
    PointIter cache) const
{
    int size = deg;
    auto t_ =  t + nu;
    assert(nu >= size_t(size));
    for(int sz = size; sz > 0; --sz) {
        for(int j = 1; j <= sz; ++j) {
            double r = t_[j] - t_[j - sz];
            double lambda =  sdiv(us[sz-1]  - t_[j - sz], r);
            cache[j - 1] = lerp(lambda, cache[j - 1], cache[j] );
        }
    }
}

// utility class derived from above that works with a vector<double> iterator
struct rmat_base_vd: rmat_base < std::vector<double>::const_iterator >
{
    typedef std::vector<double> knots_t;
    rmat_base_vd(const knots_t & t, int deg):
        rmat_base(t.cbegin(), t.cend(), deg){}

};

template<class Point>
struct rmat : public rmat_base_vd {
    typedef Point point_t;
    typedef decltype(mk_stdvec(point_t())) cpts_t;

    rmat(const cpts_t& cpts_,const knots_t& t,int deg)
        : rmat_base_vd(t, deg), cpts(cpts_)
    {
    }

    // evaluate the bspline curve at parameter u
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
    void eval_derivatives(int derOrder,
                          ParamIter us, ParamIter end,
                          PointOutIter out)
    {
        int size = degree();
        std::unique_ptr<point_t[]> cache(new point_t[size+1]);

        size_t nu = -1;
        for( ;us != end; ++us,++out) {

            nu = locate_nu( * us,nu);

            for(int j = 0; j < size + 1; ++j)
                cache[j] = control_pt(j + nu - size);

            spline_compute(nu,* us, derOrder, cache.get());
            *out =  cache[0];
        }
    }

    template <class PointOutIter>
    void eval_derivatives(int derOrder,
                          double u,
                          PointOutIter out)
    {
        int size = degree();
        std::unique_ptr<point_t[]> cache(new point_t[size+1]);
        size_t nu = locate_nu(u);

        for(int i =  0; i <= derOrder;++i,++out) {

            for(int j = 0; j < size + 1; ++j)
                cache[j] = control_pt(j + nu - size);

            spline_compute(nu, u, i, cache.get());
            *out =  cache[0];
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

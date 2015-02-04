//-*- mode:c++ -*-
#ifndef RMAT_HPP
#define RMAT_HPP
#include <vector>

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

    // find knot range [t[mu],t[mu+1]) in which @u belongs
    KnotIter
    locate(double u) const;

    double front() const {  return  t[0]; } ;
    double back() const { return   e[-1] ;};

    size_t start_mult() ;
    size_t end_mult();

    // find degree of the basis polynomials
    int
    degree() const {
        return deg;
    }

    /*const knots_t &
    knots() const {
        return t;
    }*/

    double
    coeff(size_t idx,  double u) const {
        return der_n(idx, 0, u);
    }
    double
    der(size_t idx,  double u) const {
        return der_n(idx, 1, u);
    }

    inline static double sdiv(double n, double d) {
        if (tol::param_eq(d, 0) && tol::param_eq(n, 0))
            return 0.0;
        else
            return n / d;
    }
    inline static double sdiv1(double n, double d) {
        if (tol::param_eq(d, 0) && tol::param_eq(n, 0))
            return 1.0;
        else
            return n / d;
    }

    template <class T>
    void spline_compute(size_t nu,
                        double u,
                        int numDer,
                        T cache) const
    {
        int size = deg;
        auto t_ =  t + nu;
        size_t fac = 1;

        if(numDer > deg )
            cache[0] = std::decay<decltype(cache[0])>::type(0.0);

        for(int sz = size; sz > size - numDer; --sz) {
            fac *= sz;
            for(size_t j = 1; j <= sz; ++j) {
                double d = t_[j] - t_[j - sz];
                if(d < tol::param_tol)
                    throw geom_exception(knot_not_in_range_error_der);

                double lambda =  1 / d;
                cache[j - 1] = dlerp(lambda, cache[j - 1], cache[j]);
            }
        }

        for(int sz = size - numDer; sz > 0; --sz) {
            for(size_t j = 1; j <= sz; ++j) {
                double d = t_[j] - t_[j - sz];
                double lambda =  sdiv(u  - t_[j - sz], d);
                cache[j - 1] = lerp(lambda, cache[j - 1], cache[j] );
            }
        }
        scale(cache[0],double(fac));
    }
    //.  basis_comp.png
    // returns the vector b = (B_{\mu - p, p}(u), \ldots, B_{\mu, p}(x)), mu =
    // int such that u\in [t_\mu, t_\mu + 1) and p is the degree
    std::vector<double> basis(double u) {
        int p = degree();
        size_t nu = locate_nu(u);
        std::vector<double> b(p + 1, 0.0);
        b[p] = 1;
        for(int r = 1;r <= p; ++r)
        {
            size_t k = nu - r + 1;
            double d = (t[k + r] - t[k]);
            double lambda2 = 0;
            lambda2 = sdiv(t[k + r] - u, d);
            b[p - r] = lambda2 * b[p - r + 1];
            for(int i = p - r + 1;i < p; ++i)
            {
                ++k;
                b[i] *= (1 - lambda2);
                double d = (t[k + r] - t[k]);
                lambda2 = sdiv(t[k + r] - u, d);
                b[i] += lambda2 * b[i + 1];
            }
            b[p] *= (1 - lambda2);
        }
        return b;
    }

    template <class KI, class T>
    void spline_compute(KI us,
                   size_t nu, //location of us[0] in knot seq t
                   T cache) const

    {
        int size = deg;
        auto t_ =  t + nu;
        assert(nu >= size);
        for(int sz = size; sz > 0; --sz) {
            for(int j = 1; j <= sz; ++j) {

                double r = t_[j] - t_[j - sz];

                double lambda =  sdiv(us[sz]  - t_[j - sz], r);

                cache[j - 1] = lerp(lambda, cache[j - 1], cache[j] );
            }
        }
    }

    double der_n(size_t idx,int numDer, double u) const;

    size_t locate_nu(double u) const {
        auto nu = std::distance(t, locate(u));
        return nu;
    }
    size_t
    locate_nu(double u, size_t nu_guess) const;
protected:
    KnotIter t,e;
    int deg;
};

struct rmat_base_vd: rmat_base < std::vector<double>::const_iterator >
{
    typedef std::vector<double> knots_t;
    rmat_base_vd(const knots_t & t, int deg):
        rmat_base(t.cbegin(), t.cend(), deg){}
};

template<class Point>
struct rmat : public rmat_base_vd {

    typedef Point point_t;
    typedef std::vector<point_t> cpts_t;

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
    eval_derivative(int numDer, double u) const;

    template <class ParamIter,class PointIter>
    void
    eval_pts(ParamIter f, ParamIter l, PointIter out) const {
        return eval_derivatives(0,f,l, out);
    }

    template <class ParamIter,class PointOutIter>
    void eval_derivatives(int numDer,
                          ParamIter us, ParamIter end,
                          PointOutIter out)
    {
        int size = degree();

        std::unique_ptr<point_t[]> cache(new point_t[size+1]);

        size_t nu = -1;
        for( ;us != end; ++us,++out) {

            nu = locate_nu(u,nu);

            for(int j = 0; j < size + 1; ++j)
                cache[j] = control_pt(j + nu - size);

            spline_compute(nu, u, numDer, cache.get());
            *out =  cache[0];
        }
    }

    std::vector<point_t>
    insert_knots(const std::vector<double>& taus);

    const cpts_t& control_pts() const {
        return cpts;
    }
    const point_t& control_pt(size_t i) const {
        return cpts[i];
    }
    const cpts_t& cpts;
};
}
#endif //RMAT_HPP

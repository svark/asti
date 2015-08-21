//-*- mode:c++ -*-
#include <vector>

#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/tuple.hpp>

#include <numeric>

#include "insert_knot.hpp"
#include "tol.hpp"
#include "rmat.hpp"
#include "geom_exception.hpp"
#include <random>
#include "spline_approx.hpp"
#include "bspline_cons.hpp"
#include "spline_traits.hpp"
#include "point.hpp"
#include "find_param.hpp"

namespace geom {
//{{{ -- approximations

template <class Fn>
bspline<double>
ops::cubic_approx1d(Fn f, std::vector<double> t)
{

    double mindist = std::numeric_limits<double>::infinity();
    size_t i  = 0;
    static const int p = 3;
    size_t n = t.size() - p - 1;
	assert(n>0);
    std::vector<double> pts(n,0.0);
    // pg 173 lyche
    //(@file :file-name "./lee_quasi.pdf" :to "./lee_quasi.pdf" :display "intro to quasi interpolants")

    Eigen::Matrix<double, 5, 5, Eigen::RowMajor>  mat(5, 5);
    mat.setZero();
    rmat_base_vd rm(t, 3);
    for(size_t j = 2; j + 2 < n; ++j) {
        double t5[] = { t[j + 1], (t[j + 1] + t[j + 2])/2,
                        t[j + 2], (t[j + 2] + t[j + 3])/2,
                        t[j + 3] - tol::param_tol/2 }; // ./media/cubic_approx0.png

        for(int i =  0; i < 5; ++i) {
            auto const & b = rm.basis(t5[i]).first;
            for(int k = 0; k < 4; ++k) {// ./media/cubic_approx.png
                mat(i, k + (i > 1? 1 : 0)) = b[k];
            }
        }
        double rhs[5] = {f(t5[0]), f(t5[1]), f(t5[2]), f(t5[3]), f(t5[4])};
        Eigen::Matrix<double, 5, 1> rhs_(rhs);
        Eigen::Matrix<double, 5, 1> res = mat.lu().solve(rhs_);
        if(j==2) {
            pts[j-2] = res(0);
            pts[j-1] = res(1);
        }
        pts[j] = res(2);
        if(j == n-3) {
            pts[j+1] = res(3);
            pts[j+2] = res(4);
        }
    }

    return make_bspline(std::move(pts), std::move(t), p);
}

template <class FnType>
bspline<double>
ops::quad_approx1d(FnType f, std::vector<double> t_)
{
    assert(t_.size()!=0);
    auto b = t_.cbegin();
    auto e = t_.cend();

    static const int p = 2;
    size_t n = t_.size() - p - 1;
    std::vector<double> pts(n,0);

    double s = *b;
    double t;
    pts[0] = f(s);

    for(size_t j = 1; j < n - 1; ++j) {
        s = b[j];
        t = b[j+1];
        auto u = (s + t)/2;
        pts[j] = ( -f(s) + 4*f(u) - f(t) )/2;
    }

    t = e[-1];
    pts[n-1] = f(t);

    return make_bspline(std::move(pts), std::move(t_), p);
}
//}}}

//{{{ --(@* "compute foot param of a point on a spline")
namespace kts{

template <int p>
static void build_knots_helper(const std::vector<double>& uniqts,
                               std::vector<double>& t)
{
    assert(uniqts.size() >= 2);

    std::vector<double> taus;
    taus.reserve(2*p+2);
    std::vector<int> indices(uniqts.size() - 1);
    std::iota(indices.begin(), indices.end(), 0);

    auto width_comp = [&uniqts](int i, int j) -> bool {
        double w = uniqts[i+1] - uniqts[i];
        double u = uniqts[j+1] - uniqts[j];
        return w < u;
    };

    std::make_heap(indices.begin(), indices.end(),
                   width_comp);

    std::mt19937 gen(1013); // todo:replace with linear congruence
    std::uniform_real_distribution<> rb(0.0, 1.0);
    while( uniqts.size() + taus.size() < p + 1 ) {
        int indx_ = indices[0];
        std::pop_heap(indices.begin(),
                      indices.end(),
                      width_comp);

        auto   it     = uniqts.cbegin() + indx_;
        auto   upIt   = std::next(it);
        double lambda = rb(gen);
        taus.push_back(lambda*(*upIt) + (1 - lambda)*(*it));
    }
    std::sort(taus.begin(), taus.end());
    std::merge(uniqts.begin(), uniqts.end(), taus.begin(),
               taus.end(), std::back_inserter(t));
	return ;
}

template <int p, class KnotIter>
static void build_knots(KnotIter b,
                                  KnotIter e,
                                  std::vector<double>& t,
                                  periodic_tag)
{
    size_t sz = std::distance(b, e);
    t.reserve(std::max(sz, size_t(2*p)));
    std::vector<double> uniqts;
    uniqts.reserve(sz);
    std::unique_copy(b, e, std::back_inserter(uniqts),
                     tol::param_eq );

    return build_knots_helper<p>(uniqts, t);
}

template <int p, class KnotIter>
static void build_knots(KnotIter b,
                        KnotIter e,
                        std::vector<double>& t, regular_tag)
{
    size_t sz = std::distance(b, e);
    t.reserve(sz + 2*p);

    std::vector<double> uniqts;
    uniqts.reserve(sz);
    std::unique_copy(b, e, std::back_inserter(uniqts),
                     tol::param_eq );

    for(int i = 0; i < p;++i)
        t.push_back(*b);

    build_knots_helper<p>(uniqts, t);

    for(int i = 0; i < p;++i)
        t.push_back(e[-1]);

	return ;
}
}



template <class SplineType>
double
ops::foot_param(const SplineType &spl,
                const typename SplineType::point_t& p)
{
    auto pr = spl.param_range();
    double b = pr.first;
    double e = pr.second;

    auto dist = [&p, &spl] (double u) -> double {
        return sqlen( p - spl.eval(u) );
    };

    auto dist_der_by_2 = [&p, &spl] (double u) -> double  {
        return dot( spl.eval(u) - p,  spl.eval_derivative(1, u) );
    };
    namespace m = boost::math;
    auto fn = [&p, &spl]
        (double u) {
        auto v         =  spl.eval(u) - p;
        auto vdash     =  spl.eval_derivative(1, u);
        auto vdashdash =  spl.eval_derivative(2, u);
        return m::make_tuple(dot(v, vdash),
                             dot(v, vdashdash ) + dot(vdash, vdash),
                             3*dot(vdash, vdashdash)
            );
    };

    auto deg  = spl.degree();
    auto &t = spl.knots();

    std::vector<double> tapprox;
    typedef typename spline_traits<SplineType>::ptag ptag;
    kts::build_knots<3>(spl.knots().cbegin() + spl.degree(),
                        spl.knots().cend()   - spl.degree(),
                        tapprox,
                        ptag());

	
    // get a spline approximation of the derivative of distance
    auto splapprox = cubic_approx1d(dist_der_by_2, tapprox);

    // compute roots of this spline
    auto rootQ = find_next_rootc(splapprox, t.front(), tol::param_tol);

    double mindist = dist(t.front());
    double minu = t.front();
    double min = minu;

    for(;rootQ.second;) {
        const size_t digits_ = 2 *std::numeric_limits<double>::digits/3;

        double guess = rootQ.first;

        double u = m::tools::halley_iterate(fn, guess,
                                            minu + tol::param_tol,
                                            t.back(), //TODO:adjust
                                            digits_);
        min = std::max(min, u);
        auto d = dist(u);
        if( d < mindist ) { mindist = d; minu = u; }
        rootQ = find_next_rootc(splapprox,
                                std::max(u, rootQ.first)+tol::param_tol,
                                tol::param_tol);//TODO:adjust
    }
    double db = dist(t.back());
    if(mindist > db ) {
        minu = t.back();
    }
    return minu;
}
//}}}
}

/*
  Local Variables:
  eval:(load-file "./scripts/temp.el")
  eval:(setq methods (list "foot_param"
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
  eval:(instantiate-templates "spline_approx" "ops" (list ) (product methods spltypes) )
  End:
*/
#include "bspline.hpp"
#include "periodic_bspline.hpp"
#include "rational_bspline.hpp"
#include "point.hpp"
#include <functional>
namespace geom {
#include "spline_approx_inst.inl"
template bspline<double> ops::cubic_approx1d(
    class std::function<double (double)>, class std::vector<double> );
template bspline<double> ops::quad_approx1d(
    class std::function<double (double)>, class std::vector<double> );
}

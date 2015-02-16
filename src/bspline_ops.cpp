#include <stdafx.h>

#include <random>
#include <functional>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>

#include "boost/mpl/integral_c.hpp"
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/tuple.hpp>
#include "geom_exception.hpp"
#include "bspline.hpp"
#include "periodic_spline.hpp"
#include "bspline_ops.hpp"
#include "rmat.hpp"
#include "tol.hpp"
#include "util.hpp"
#include "spline_traits.hpp"
#include "smat.hpp"
namespace geom {


//{{{ --(@* "implements the oslo algorithm to insert knots")
// knots @us each of which repeated as in @numRepeats
// see(@url :file-name "http://www.uio.no/studier/emner/matnat/ifi/INF-MAT5340/v09/undervisningsmateriale/book.pdf#page=100" :display "page100")
template <class SplineType>
static SplineType
bspline_ops::insert_knots(const SplineType& crv,
                          const std::vector<double>& new_knots)
{
    typedef typename SplineType::point_t point_t;
    typedef typename SplineType::knots_t knots_t;
    auto const &uarr = new_knots;

    knots_t taus;
    auto const & t = crv.knots();

    taus.reserve( uarr.size() + t.size() );

    auto const & cpts = crv.control_points();
    std::merge(t.begin(), t.end(),
               uarr.begin(), uarr.end(),
               std::back_inserter(taus));

    rmat<point_t> m(cpts, t, crv.degree());
    auto new_cpts(m.insert_knots(taus));
    return SplineType(std::move(new_cpts),
                      std::move(taus),
                      crv.degree()).translate(crv.base_point());
}

template <class SplineType>
static SplineType
bspline_ops::insert_knot(const SplineType& crv,
                         double u)
{
    typedef typename SplineType::point_t point_t;
    typedef typename SplineType::knots_t knots_t;

    double uarr[] = {u};
    auto const & t = crv.knots();
    knots_t taus( 1 + t.size() );
    taus.assign(t.begin(), t.end());
    taus.insert(std::upper_bound(taus.begin(), taus.end(), u), u );
    auto const & cpts = crv.control_points();

    rmat<point_t> m(cpts, t, crv.degree());
    auto new_cpts(m.insert_knots(taus));
    return SplineType(std::move(new_cpts),
                      std::move(taus),
                      crv.degree()).translate(crv.base_point());
}


//}}}

//{{{ --(@* "raise the degree by 1")
// see:(@file :file-name "raise_degree.png" :to "./raise_degree.png" :display "raise_degree")
template <class SplineType>
SplineType bspline_ops::raise_degree(const SplineType& spl)
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
        auto n = f;
        // increase the mulitplicity of each knot by 1
        if(++n != e && !tol::param_eq(*n, *f))
            new_knots.push_back(*f);
    }

    size_t num_new_knots = new_knots.size();
    int p = degree() + 1;
    size_t num_new_cpts = num_new_knots - p - 1;

    using util::skip_ith_iter;

    for(size_t i = 0; i < num_new_cpts; ++i)
    {
        vector_t cv(0);
        for(size_t j = 0; j < p; ++j) {
            skip_ith_iter iter( j, new_knots.cbegin() + i );
            cv += make_vec(spl.blossom_eval(iter)) ;
        }
        cv *= 1.0/p;
        new_cpts.push_back(make_pt(v));
    }
    return SplineType(std::move(new_cpts),
                      std::move(new_knots), p).translate(spl.base_point());
}
//}}}

//{{{ --(@* "reverse curve sense")
template <class SplineType>
static SplineType bspline_ops::reverse_curve(SplineType& spl)
{
    auto &cpts=spl.control_points();
    auto &t = spl.knots();
    size_t numPts = cpts.size();
    typedef typename SplineType::point_t point_t;
    std::vector<point_t> new_cpts(numPts);
    new_cpts.reserve(numPts);

    std::vector<double> new_knots(spl.knots().size());

    std::reverse_copy(cpts.cbegin(), cpts.cend(), new_cpts.begin());
    std::reverse_copy(t.cbegin(), t.cend(), new_knots.begin());

    std::transform(new_knots.begin(), new_knots.end(),
                   new_knots.begin(),
                   std::negate<double>() );

    return SplineType(std::move(new_cpts), std::move(new_knots),
                      spl.degree()).translate(spl.base_point());
}

template <class SplineType>
static SplineType& bspline_ops::inplace_reverse_curve(SplineType& spl)
{
    //typedef typename SplineType::point_t point_t;
    std::reverse(spl.control_points().begin(), spl.control_points().end());
    std::reverse(spl.knots().begin(), spl.knots().end());

    std::transform(spl.knots().begin(), spl.knots()().end(),
                   spl.knots().begin(),
                   std::negate());
    return spl;
}

//}}}

//{{{ --(@* "re-parametrise curve with a new range")
template <class SplineType>
static SplineType bspline_ops::reparametrize(const SplineType& spl,
                                             double t1, double t2)
{
    typedef typename SplineType::knots_t knots_t;
    typedef typename SplineType::point_t point_t;
    knots_t new_knots( spl.knots().size() );
    auto & t = spl.knots();

    double last_t = t.back();
    double first_t = t.front();

    if( fabs(last_t - first_t) < tol::param_tol)
        throw geom_exception(bad_knot_spacing_t);

    const double scale = 1.0 / ( last_t - first_t);
    size_t i = 0;
    for(double u : t )
    {
        const double par = ( u - first_t) * scale;
        new_knots[i++] = t1 + (t2 - t1) * par;
    }
    return SplineType(
          SplineType::cpts_t(spl.control_points()),
          std::move(new_knots),
          spl.degree()
        ).translate(spl.base_point());
}


//}}}

//{{{ --(@* "is bezier")
template <class SplineType>
static bool bspline_ops::is_bezier(const SplineType& spl)
{
    auto &t = spl.knots();
    int sz = spl.degree() + 1;
    if( t.size() != 2*sz )
        return false;
    if(std::equal(t.begin(), t.begin() + sz,
                   t.rbegin() + t.size() - sz,  tol::param_eq ))
    {
        if(std::equal(t.rbegin(), t.rbegin() + sz,
                   t.begin() + t.size() - sz,  tol::param_eq ))
        {
            return true;
        }
    }
    return false;
}
//}}}

//{{{ --(@* "split into bezier patches")

template <class SplineType>
static std::list<SplineType>
bspline_ops::split_into_bezier_patches(const SplineType &spl)
{
    typedef typename SplineType::knots_t knots_t;
    typedef typename SplineType::point_t point_t;

    auto &ts    = spl.knots();
    int   deg   = spl.degree();
    int   order = deg + 1;

    assert(!tol::param_eq(ts.front(), ts.back()));
    SplineType s ( clamp_end(clamp_start(spl)) );
    std::list<SplineType> patches;
    knots_t newts;
    newts.reserve( s.knots().size() );
    std::unique_copy(s.knots().begin(),
                     s.knots().end(),
                     std::back_inserter(newts),
                     tol::param_eq);

    for(int j = 1;j < newts.size()-1; ++j)
    {
        patches.push_back(clamp_at_right(newts[j], s));
        s.swap(clamp_at_left(newts[j], s));
    }
    patches.push_back(s);
	return patches;
}

template <class SplineType>
static SplineType
bspline_ops::extract_curve(const SplineType &spl, double a, double b)
{
    typedef typename SplineType::knots_t knots_t;
    typedef typename SplineType::point_t point_t;

    auto &ts    = spl.knots();
    int   deg   = spl.degree();
    int   order = deg + 1;

    assert(!tol::param_eq(ts.front(), ts.back()));
    return clamp_end(clamp_start(spl));
}

template <class SplineType>
static std::list<SplineType>
bspline_ops::split_into_bezier_patches_hard(const SplineType &spl)
{
    typedef typename SplineType::knots_t knots_t;
    typedef typename SplineType::point_t point_t;

    auto &ts    = spl.knots();
    int   deg   = spl.degree();
    int   order = deg + 1;

    assert(!tol::param_eq(ts.front(), ts.back()));

    knots_t newts;
    newts.reserve( order * ts.size() );
    // newts.push_back(ts.front());

    for(auto it = std::next(ts.begin()); it != std::prev(ts.end()); ++it)
    {
        int mult = 1;
        double u = *it;
        while(std::next(it)!= std::prev(ts.end()) &&* std::next(it) == *it)
        {
            ++mult;
            ++it;
        }
        for(;mult != deg; ++mult)
            newts.push_back(u);
    }

    // newts.push_back(ts.back());
    auto refined_spl(bspline_ops::insert_knots(spl, newts));
    auto const &cpts = refined_spl.control_points();

    std::list<SplineType> patches;

    //  {0, 0, 0, 1    , 2    , 3, 3, 3} -->
    //  {0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3}

    //  patch1 ts- (0, 0, 0, 1, 1, 1),
    //  with cpts B(0, 0, 0), B(0, 0, 1), B(0, 1, 1), B(1, 1, 1)

    //  patch2 ts- (1, 1, 1, 2, 2, 2)
    //  with cpts B(1, 1, 1), B(1, 1, 2), B(1, 2, 2), B(2, 2, 2)

    for(size_t sz = 0;sz < cpts.size(); sz += deg ) {

        std::vector<point_t>  patch_cpts( cpts.begin() + sz,
                                          cpts.begin() + deg + 1 + sz);

        knots_t  patch_ts ( refined_spl.knots().cbegin() + sz,
                            refined_spl.knots().cbegin() + sz + 2 * deg);

        patches.push_back( SplineType( std::move(patch_cpts),
                                       std::move(patch_ts), deg)
                           .translate(spl.base_point()) );
    }
    return patches;
}


//}}}

//{{{ --
template <class PointIter, class KnotIter>
bool bspline_ops::is_periodic( PointIter pb, PointIter pe,
                               KnotIter tb, KnotIter te,
                               int deg)
{
    int p = deg;
    size_t np = std::distance(pb, pe);
    size_t nt = std::distance(tb, te);
    if( np < 2)
        return false;

    if( nt != np + p)
        return false;

    if(nt < 2*p + 2 )
        return false;

    if( pb[0] != pe[-1] )
        return false;

    std::vector<double> buf1(p + 1), buf2(p + 1);

    te -= (p + 1);

    auto deltas_at_start = buf1.begin();
    // fill up (1 + p) knot ranges) at the start into buf1
    std::adjacent_difference(tb,
                             tb + p + 1,
                             deltas_at_start);

    if( std::all_of( deltas_at_start + 1,
                     deltas_at_start + p + 1,
                     [](double v)->bool{ return v == 0; }  ))
        return false;

    auto deltas_at_end = buf2.begin();
    // fill up 1 + (p knot ranges) at the end into buf1
    std::adjacent_difference(te,
                             te + p + 1,
                             deltas_at_end);

    if( std::all_of( deltas_at_end + 1,
                     deltas_at_end + p + 1,
                     [](double v)->bool{ return v == 0; }  ))
        return false;

    return std::equal(deltas_at_start + 1,
                      deltas_at_start + p,
                      deltas_at_end + 1);
}
//}}}

//{{{ --
template <class Fn>
static bspline<double>
bspline_ops::cubic_approx1d(Fn f, std::vector<double>& t)
{

    double mindist = std::numeric_limits<double>::infinity();
    size_t i  = 0;
    static const int p = 3;
    size_t n = t.size() - p - 1;
    std::vector<double> pts(n,0.0);
    // pg 173 lyche
    //(@file :file-name "./lee_quasi.pdf" :to "./lee_quasi.pdf" :display "intro to quasi interpolants")

    Eigen::Matrix<double, 5, 5, Eigen::RowMajor>  mat(5, 5);
    mat.setZero();
    rmat_base_vd rm(t, 3);
    for(size_t j = 2; j + 2 < n; ++j) {
        double t5[] = { t[j + 1], (t[j + 1] + t[j + 2])/2,
                        t[j + 2], (t[j + 2] + t[j + 3])/2,
                        t[j + 3] - tol::param_tol/2 }; // ./cubic_approx0.png

        for(int i =  0; i < 5; ++i) {
            auto const & b = rm.basis(t5[i]);
            for(int k = 0; k < 4; ++k) {// ./cubic_approx.png
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

    return bspline<double>(std::move(pts), std::move(t), 3);
}

template <class FnType>
static bspline<double>
bspline_ops::quad_approx1d(FnType f, std::vector<double>& t)
{
    assert( std::distance(b, e)!=0);

    if(e[-1] - b < param_tol)
        throw geom_exception(bad_knot_spacing_t);

    std::vector<double> pts(n);

    double s = *b;
    double t;
    pts[0] = f(s);

    for(size_t j = 1; j < n - 1; ++j) {
        s = b[j];
        t = b[j+1];
        u = (s + t)/2;
        pts[j] = ( -f(s) + 4*f(u) - f(t) )/2;
    }

    t = e[-1];
    pts[n-1] = f(t);

    return bspline<double>(std::move(pts), std::move(knots), 3);
}
//}}}

//{{{ --(@* "compute foot param of a point on a spline")
namespace kts{

template <int p>
void build_knots_helper(const std::vector<double>& uniqts,
                        std::vector<double>& t)
{
    if(uniqts.size() < 2)
        throw geom_exception(bad_knot_spacing_t);

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
}

template <int p, class KnotIter>
void build_knots(KnotIter b,
                 KnotIter e,
                 std::vector<double>& t, periodic_tag)
{
    size_t sz = std::distance(b, e);
    t.reserve(std::max(sz, size_t(2*p)));
    std::vector<double> uniqts;
    uniqts.reserve(sz);
    std::unique_copy(b, e, std::back_inserter(uniqts),
                     tol::param_eq );

    build_knots_helper<p>(uniqts, t);
}

template <int p, class KnotIter>
void build_knots(KnotIter b,
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
}
}

std::pair<double, bool>
find_next_rootc(bspline<double>& spl,
                double prev,
                double tol)
{
    auto const &t = spl.knots();
    auto const &c = spl.control_points();
    int p = spl.degree();

    size_t k = 1;
    while( k < c.size() && t[k] < prev)
        ++k;

    for(; ; )
    {
        auto const &t = spl.knots();
        auto const &c = spl.control_points();
        auto tcap = [&t, p](size_t i) -> double {
            return std::accumulate( &t[i+1], &t[i + 1] + p, 0.0 )/p;
        };

        size_t n = c.size();
        while(k < n && ( c[k-1]*c[k] > 0 ))
            ++k;

        if( k >= n)
            break;

        double root =
            tcap(k) - (c[k]/p) * (t[k+p] - t[k] )/(c[k] - c[k-1] );

        if( fabs(spl.eval(root)) < tol)
            return std::make_pair(root, true);

        spl.swap(bspline_ops::insert_knot(spl, root));
    }
    return std::make_pair(prev, false);
}


template <class SplineType>
static double
bspline_ops::foot_param(const SplineType &spl,
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
        auto v =  spl.eval(u) - p;
        auto vdash =  spl.eval_derivative(1, u);
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
                                std::max(u, rootQ.first)+2.0e-3,
                                tol::param_tol);//TODO:adjust
    }
    double db = dist(t.back());
    if(mindist > db ) {
        minu = t.back();
    }
    return minu;
}
//}}}

//--compute box
template <class SplineType>
static box<spline_traits<SplineType>::dim>
bspline_ops::compute_box(const SplineType &spl)
{
    box < spline_traits < SplineType >::dim > b;
    auto & cpts = spl.control_points();
    for(auto c : cpts)
    {
        b += c;
    }
    b.lo += spl.base_point();
    b.hi += spl.base_point();
    return b;
}
}

#include "point.hpp"


/*
Local Variables:
eval:(load-file "./temp.el")
eval:(setq methods (list "insert_knots"
                         "reparametrize"
                         "split_into_bezier_patches_hard"
                         "split_into_bezier_patches"
                         "extract_curve"
                         "is_bezier"
                         "reverse_curve"
                         "foot_param"
                         "compute_box"
                         ))
eval:(setq spltypes (list "bspline<double>"
                           "bspline<point2d_t>"
                           "bspline<point3d_t>"
                           "bspline<point4d_t>" ))
eval:(instantiate-templates "bspline_ops" (list ) methods spltypes )
End:
// dump all explicitly instantiated templates below
*/
namespace geom {
#include "bspline_ops_inst.cpp"
template bspline<double> bspline_ops::cubic_approx1d(class std::function<double (double)>,class std::vector<double> &);
}

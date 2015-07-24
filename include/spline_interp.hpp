#ifndef ASTI_SPLINE_INTERP_HPP
#define ASTI_SPLINE_INTERP_HPP
#include "util.hpp"
#include "Eigen/Core"
#include <Eigen/Dense>

#include "point_dim.hpp"
#include "bspline.hpp"
#include "periodic_bspline_cons.hpp"
#include <type_traits>
#include <numeric>
#include "geom_exception.hpp"
#include "tol.hpp"
#include "match_case.hpp"
namespace geom {


// {{{ -- interpolation (hermite & global)
template <int dim,class PointIter>
void
setupQMatrix(PointIter pb, PointIter pe,
             Eigen::Matrix<double,dim,
             dim >& sigmaxy)
{
    typedef RAWTYPE(pb[0]) point_t;
    point_t cg(centroid(pb, pe));
    sigmaxy.setZero();
    //should be able to run in parallel
    int n = 0;
    for(;pb!=pe;++pb)
    {
        vec_t<dim> p(*pb-cg);
        for(int i = 0;i < dim;++i)
        {
            for(int j = i ; j < dim;++j)
                sigmaxy(i,j) += (p[i])*(p[j]);
        }
        ++n;
    }

    for(int i = 0; i < dim;++i) {
        for(int j = i; j < n ; ++j){
            sigmaxy(i,j) /= n;
            sigmaxy(j,i) = sigmaxy(i,j);
        }
    }
}

enum parametrization_option_t { centripetal_length, chord_length, affinely_invariant,neilson_foley};
enum end_conditions_t { periodic, parabolic_blending,vanishing_double_derivatives,not_a_knot} ;
struct interpolation_options_t
{
    parametrization_option_t  parametrization;
    end_conditions_t  end_conditions;
};


template <class PointIter,class KnotIter>
void find_parameters( PointIter pb, PointIter pe,
                      KnotIter tb,
                      std::integral_constant<parametrization_option_t,
                      chord_length> )
{
    std::size_t n = std::distance(pb,pe);
    std::vector<double> widths(n);
    widths[0] = 0.0;
    for (size_t i = 0; i < n - 1; ++i) {
        double dist = 0, l = 0;
        auto v = pb[i + 1] - pb[i];
        l = sqrt(dot(v,v));
        widths[i+1] = l;
    }
    std::partial_sum(widths.begin(),widths.end(),tb);
}

template <class PointIter,class KnotIter>
void find_parameters( PointIter pb, PointIter pe,
                      KnotIter tb,
                      std::integral_constant<parametrization_option_t,
                      centripetal_length> )
{
    std::size_t n = std::distance(pb,pe);
    std::vector<double> widths(n);
    widths[0] = 0.0;
    for (size_t i = 0; i < n - 1; ++i) {
        double l = 0;
        auto v = pb[i + 1] - pb[i];
        l = sqrt(sqrt(dot(v,v)));
        widths[i+1] = l;
    }
    std::partial_sum(widths.begin(),widths.end(),tb);
}

template <class PointIter,class KnotIter>
void find_parameters(PointIter pb,
                     PointIter pe,
                     KnotIter t,
                     std::integral_constant<parametrization_option_t,
                     affinely_invariant>)
{
    std::size_t n = std::distance(pb,pe);
    std::vector<double> widths(n);
    widths[0] = 0.0;
    static const int dim = point_dim<RAWTYPE(pb[0])>::dimension;
    Eigen::Matrix<double,dim,dim> Q;
    setupQMatrix<dim>(pb,pe,Q);
    for (size_t i = 0; i < n - 1; ++i) {
        double l = 0;
        vec_t<dim> diff(  pb[i + 1] - pb[i] );
        auto  v = eigen_vec ( diff );
        l = sqrt (v.transpose() * (Q) *  v) ;
        widths[i+1] = l;
    }
    std::partial_sum(widths.begin(),widths.end(),t);
}

template <class PointIter,class KnotIter>
void find_parameters(PointIter pb, PointIter pe,
                     KnotIter t,
                     std::integral_constant<parametrization_option_t,
                     neilson_foley> )
{
    std::size_t n = std::distance(pb,pe);
    std::vector<double> widths(n), thetas(n);
    widths[0] = 0.0;
    thetas[0] = 0.0;
    static const int dim = point_dim<RAWTYPE(pb[0])>::dimension;
    Eigen::Matrix<double,dim,dim> Q;
    setupQMatrix<dim>(pb,pe,Q);
    RAWTYPE(eigen_vec(pb[0])) u;
    for (size_t i = 0; i < n - 1; ++i) {
        double l = 0;
        vec_t<dim> diff(  pb[i + 1] - pb[i] );
        auto const & v = eigen_vec(diff);
        l = sqrt( (v.transpose() * (Q) *  v)(0) );
        widths[i+1] = l;
        if(i > 0 )
            thetas[i] = acos( (u.transpose()*v)(0) /(u.norm()*v.norm()) );
        u = v;
    }

    for (size_t i = 1; i < n - 1; ++i)
        widths[i] *= (1
                     + 3 *thetas[i] * widths[i - 1]
                     / 2*( widths[i - 1] + widths[i] )
                     + 3 *thetas[i + 1] * widths[i + 1]
                     / 2*( widths[i] + widths[i + 1] ) );

    std::partial_sum(widths.begin(),widths.end(),t);
}

template <class Matrix,class PointIter,
          class KnotIter,
          class Tangents, class MatrixOut>
void setup_mat_for_interior_tgts(Matrix &m,
                                 PointIter  p, PointIter  pe,
                                 KnotIter   t, KnotIter   te,
                                 const Tangents& vecs,
                                 MatrixOut& d )
{
// set up rows 1 through n-2 of the matrix
// as in hoschek  pg 88
// eqn 3.18. contd.
    size_t n = std::distance(t,te);
    assert( n > 2);
    static const int dim = point_dim<RAWTYPE(p[0])>::dimension;

    for (size_t j = 1; j < n - 1; ++j) {
        auto aj = E(t,j);
        auto aj1 = E(t,j-1);

        m(j, j - 1) = aj;
        m(j, j)     = 2 * (aj + aj1);
        m(j, j + 1) = aj1;

        double a = -3 * aj / aj1;
        double b =  3 * (aj / aj1  - aj1 / aj);
        double c =  3 * (aj1 / aj);

        d(j) = ( a *(p[j - 1]) + b*(p[j]) +
            c * (p[j + 1]));
    }

    // use wherever available, explicit information about derivatives.
    for(int j = 0; j < vecs.rows() ; ++j)
    {
        auto v = vecs.row(j);
        if( !tol::eq(v.norm(), 0) )
        {
            m.row(j).setZero();
            m(j, j) = 1;
            d.row(j) = (v);
        }
    }
}


template <class PointIter,class KnotIter,class VecT>
void
eval_tangents_for_pchip( PointIter pb, PointIter pe,
                         KnotIter tb, KnotIter te,
                         VecT& explicit_tgts,
                         std::integral_constant<end_conditions_t,periodic>
    )
{
    static const int dim = point_dim<RAWTYPE(pb[0])>::dimension;
    using Eigen::MatrixXd;
    size_t n = std::distance(pb,pe);

    if(!tol::eq(sqlen(pb[0]-pb[n-1]),0,1e-16))
        throw geom_exception(interp_ends_dont_match_for_periodic);

    MatrixXd mat(n,n);
    mat.setZero();
    Eigen::Matrix<double,Eigen::Dynamic,dim> rhs(n);

    mat(0,0)   = E(tb,n-3);
    mat(0,n-3) = E(tb,n-2);
    mat(0,n-2) = 2*(mat(0,0)  + mat(0,n-3));

    double a0 = (mat(0,n-3)/mat(0,0));
 
    rhs.row(0)(0) = -3*a0 * (E(pb, 0, n-2))
                 +3/a0 * (E(pb, n-3));

    mat(n-1,n-2) = E(tb,0);
    mat(n-1,1)   = E(tb,n-2);
    mat(n-1,0)   = 2*( mat(n-1,n-2) + mat(n-1,1) );

    a0 = mat(n-1,n-2)/mat(n-1,1);
    rhs.row(n-1)(0) = 3*a0* E(pb,0) - 3/a0 * E(pb,0,n-1);

    setup_mat_for_interior_tgts(mat,pb,pe,tb,te,explicit_tgts,rhs);

    rhs = mat.fullPivLu().solve(rhs);

    for(size_t i = 0; i < n; ++i)
        for(int j = 0; j < dim;++j)
            explicit_tgts.row(i).col(j) = rhs.row(i).col(j);
}

template <class PointIter,class KnotIter,class VecT>
void
eval_tangents_for_pchip( PointIter pb, PointIter pe,
                         KnotIter tb, KnotIter te,
                         VecT& explicit_tgts,
                         std::integral_constant<end_conditions_t,parabolic_blending> )
{
    static const int dim = point_dim<RAWTYPE(pb[0])>::dimension;
    size_t n = std::distance(pb,pe);
    Eigen::MatrixXd m(n,n);
    m.setZero();
    Eigen::Matrix<double,Eigen::Dynamic,dim> rhs(n);
    m(0,0)     = 1;
    m(n-1,n-1) = 1;
    //see 3.19 expanded hoschek
    rhs.row(0)   = eigen_vec((-3* pb[0] + 4*pb[1] - pb[2])/2);
    rhs.row(n-1) = eigen_vec((pb[n-3] - 4*pb[n-2] + 3*pb[n-1])/2);

    setup_mat_for_interior_tgts(m, pb, pe, tb, te, explicit_tgts, rhs);
    rhs = m.fullPivLu().solve(rhs);

    for(size_t i = 0; i < n; ++i)
        for(int j = 0; j < dim;++j)
            explicit_tgts.row(i).col(j) = rhs.row(i).col(j);
}

template <class PointIter,class KnotIter,class VecT>
void
eval_tangents_for_pchip( PointIter pb, PointIter pe,
                         KnotIter tb, KnotIter te,
                         VecT& explicit_tgts,
                         std::integral_constant<end_conditions_t,vanishing_double_derivatives> )
{
    size_t n = std::distance(pb,pe);
    Eigen::MatrixXd m(n,n);
    m.setZero();
    static const int dim = point_dim<RAWTYPE(pb[0])>::dimension;
    Eigen::Matrix<double,Eigen::Dynamic,dim> rhs(n);
    // pg 88 of hoschek..  set up equation m p' = n.p for C^2
    // continuity of cubic polynomials at params. eq 3.16 we need to
    // solve for the tangent p' ( n of them) we have assumed that at
    // param 0 and 1 the double derivatives vanish (no curvature
    // condition)
    double a0 = E(pb,0);
    double a1 = E(pb,1);
    m(0, 0) = 2 * a0;
    m(0, 1) = a0;
    a1 = E(tb,n-2);
    m(n - 1, n - 2) = a1;
    m(n - 1, n - 1) = 2 * a1;     // there is a mistake here in hoschek

    rhs.row(0)     = eigen_vec(-3 * pb[0] + 3 * pb[1]);
    rhs.row(n - 1) = eigen_vec(-3 * pb[n - 2] + 3 * pb[n - 1]);
    setup_mat_for_interior_tgts(m, pb,
                                pe, tb, te,
                                explicit_tgts, rhs);
    rhs = m.fullPivLu().solve(rhs);

   for(size_t i = 0; i < n; ++i)
        for(int j = 0; j < dim;++j)
            explicit_tgts.row(i).col(j) = rhs.row(i).col(j);
}

template <class PointIter,class KnotIter,class VecT>
void
eval_tangents_for_pchip( PointIter pb, PointIter pe,
                         KnotIter tb, KnotIter te,
                         VecT& explicit_tgts,
                         std::integral_constant<end_conditions_t,not_a_knot> )
{
    static const int dim = point_dim<RAWTYPE(pb[0])>::dimension;
    size_t n = std::distance(pb,pe);
    Eigen::MatrixXd m(n,n);
    m.setZero();
    Eigen::Matrix<double,Eigen::Dynamic,dim> rhs(n);
    auto a0 = E(tb,0);
    auto a1 = E(tb,1);
    m(0,0) = (1/a0)* ( 1/a0  +1/a1) ;

    m(0,1) = (1/a0 + 1/a1);// square this
    m(0,1) *= m(0,1);

    auto b0 = 2/(a0 * a0 * a0);
    auto b1 = 3/(a1 * a0 * a0);
    auto b2 = 1/(a1 * a1 * a1);

    rhs.row(0) = eigen_vec(
        -(b0 + b1)        * pb[0]
        +(b0  - b2  + b1) * pb[1]
        + b2              * pb[2] );

    auto h0 = E(tb,n-3);
    auto h1 = E(tb,n-2);

    m(n-1,n-2) = (h0+h1)*(h0+h1)/h1;
    m(n-1,n-1) = (h0 + h0*h0/h1);

    b0 = h0/h1;
    b1 = h1/h0;

    rhs.row(n-1) = eigen_vec( 
        - b1                  * pb[n-3]
        - (b0 + b1 + 2*b0*b0) * pb[n-2]
        + (2*b0*b0 + 3*b0   ) * pb[n-1] );

    setup_mat_for_interior_tgts(m, pb, pe, tb, te, explicit_tgts, rhs);
    rhs = m.fullPivLu().solve(rhs);

   for(size_t i = 0; i < n; ++i)
        for(int j = 0; j < dim;++j)
            explicit_tgts.row(i).col(j) = rhs.row(i).col(j);
}

void interp_fail()
{
    throw geom_exception(unknown_interp_option_t);
}

template <class PointIter, class VecsT,class ParamsT>
auto
pchip_open(PointIter pb, PointIter pe,
           const ParamsT& params, const VecsT& tgts)
           -> bspline<RAWTYPE(pb[0])>
{
    static const int dim = point_dim<RAWTYPE(pb[0])>::dimension;
    typedef RAWTYPE(pb[0]) point_t;
    size_t n = std::distance(pb,pe);
    RAWTYPE(mk_stdvec(pb[0])) cpts(2*n);
    std::vector<double> knots( 2* n + 4);

    cpts[0] = pb[0];

    cpts[1] = pb[0] + 1/3 * tgts[0] * E(params,0);

    knots[0] = knots[1] = knots[2] = knots[3] = params[0];

    for( size_t i = 2; i < n ; ++i) {
        cpts[2*i-2] = pb[i-1] - 1/3 * tgts[i-1] *  E(params,(i-2)) ;

        cpts[2*i-1] = pb[i-1] + 1/3 * tgts[i-1] *  E(params,(i-1)) ;

        knots[2*i] = knots[2*i+1] = params[i-1];
    }

    cpts[2*n-2] = pb[n-1] - 1/3 * tgts[n-1] * E(params,n-2);
    cpts[2*n-1] = pb[n-1];
    knots[2*n+3] = knots[2*n+2] = knots[2*n+1] = knots[2*n] = params[n-1];

    return make_bspline( std::move(cpts), std::move(knots ), 3);
}

template <class PointIter, class VecsT,class ParamsT>
auto
pchip_closed(PointIter pb, PointIter pe,
             const ParamsT& params,
             const VecsT& tgts
    ) -> periodic_bspline<RAWTYPE(pb[0]) >
{
    static const int dim = point_dim<RAWTYPE(pb[0])>::dimension;
    typedef RAWTYPE(pb[0]) point_t;
    size_t n = std::distance(pb,pe);
    RAWTYPE(mk_stdvec(pb[0])) cpts(2*n-1);//last cpt is omitted it is the same as the first
    std::vector<double> knots( 2* n);
//pg 180,188 hoschek
    cpts[0] = pb[0] - 1/3 * tgts[0] *  E(params,n-2) ;

    cpts[1] = pb[0] + 1/3 * tgts[0] *  E(params, 0) ;

    knots[0] = knots[1] = params[0];

    for( size_t i = 2; i < n ; ++i ) {
        cpts[ (2*i-2) ] = pb[i-1]- 1/3 *tgts[i-1] * E(params,(i-2));

        cpts[ (2*i-1) ] = pb[i-1] + 1/3 * tgts[i-1] * E(params,(i-1)) ;

        knots[2*i-2] = knots[2*i-1] = params[i-1];
    }

    cpts[2*n-2] = pb[n-1] - 1/3 * tgts[n-1] * E(params,n-2) ;

    auto cp = cpts[0] - pb[n-1] + 1/3 * tgts[n-1] * E(params,0);
    assert(tol::eq(sqlen(cp),0));

    knots[2*n-2] = knots[2*n-1] = params[n-1];

    assert( pb[n-1] == pb[0]);
    return  make_periodic_bspline_wrap(cpts, knots, 3);
}


template <class PointIter, class VecsT>
auto piecewise_cubic_hermite_interp_periodic
(
    PointIter pb,
    PointIter pe,
    interpolation_options_t opts,
    const VecsT & vecs
    ) -> periodic_bspline< RAWTYPE(pb[0]) >
{
    typedef RAWTYPE(pb[0]) point_t;
    static const int dim =  point_dim<point_t>::dimension;

    // generate parameters using chord length approximation.
    const size_t m = std::distance(pb,pe);
    assert( opts.end_conditions !=
                  periodic
                  || pe[-1] == pb[0] );
    std::vector < double > params(m);

    switch_case_4(parametrization_option_t,
                  opts.parametrization,
                  centripetal_length,
                  interp_fail,
                  find_parameters,
                  pb, pe, params.begin());

    double sum = params.cend()[-1];
    std::transform(params.cbegin(),params.cend(),params.begin(),
        [&sum](double v){ return v/sum;} );

    typedef Eigen::Matrix<double,Eigen::Dynamic,dim> MatrixXdim;

    MatrixXdim tgts(m, dim);
    tgts.setZero();

// end conditions depend on options passed..
    eval_tangents_for_pchip(
        pb,pe,params.begin(),params.end(),tgts,
        std::integral_constant<end_conditions_t,periodic>());

    return pchip_closed(pb,pe,params,tgts);
}

template <class PointIter, class VecsT>
auto
piecewise_cubic_hermite_interp_regular
(
    PointIter pb,
    PointIter pe,
    interpolation_options_t opts,
    const VecsT & vecs
    ) -> bspline< typename RAWTYPE(pb[0] ) >
{
    typedef RAWTYPE(pb[0]) point_t;
    static const int dim =  point_dim < point_t >::dimension;

    // generate parameters using chord length approximation.
    const size_t m = std::distance(pb,pe);
    assert( opts.end_conditions !=
                  periodic
                  || pe[-1] == pb[0] );

    std::vector < double > params(m);

    switch_case_4(parametrization_option_t,
                  opts.parametrization,
                  centripetal_length,
                  interp_fail,
                  find_parameters,
                  pb, pe, params.begin());

    double sum = params.cend()[-1];
    std::transform(params.cbegin(),params.cend(),params.begin(),
        [&sum](double v){ return v/sum;} );

    typedef Eigen::Matrix<double,Eigen::Dynamic,dim> MatrixXdim;

    MatrixXdim tgts(m, dim);
    tgts.setZero();

// end conditions depend on options passed..
    switch_case_5( end_conditions_t,
        opts.end_conditions, periodic, interp_fail,
        eval_tangents_for_pchip,
        pb,pe,params.begin(),params.end(),tgts);

    bool is_periodic = opts.end_conditions == periodic;
    return pchip_open(pb,pe,params,tgts);

}


template <class PointIter, class VecsT>
auto
piecewise_cubic_hermite_interp
(
    PointIter pb,
    PointIter pe,
    interpolation_options_t opts,
    const VecsT & vecs
    ) -> bspline< RAWTYPE(pb[0] ) >

{
    if(opts.end_conditions != periodic)
        return piecewise_cubic_hermite_interp_regular(pb, pe, opts, vecs);
    else
        return piecewise_cubic_hermite_interp_periodic(pb, pe, opts, vecs).spline();
}
// }}}
// }}}
}
#endif

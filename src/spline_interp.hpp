#include "util.hpp"
#include "Eigen/Core"
#include <boost/range/iterator_range_core.hpp>
namespace geom {

template < class cpts_t, class param_t  >
decltype( pts[0] - pts[1] )
start_vec(const cpts_t & pts,
          const param_t& t, int deg)
{
    typdef decltype( pts[0] - pts[1] ) vec_t;
    auto dt0 = E(t,0,deg);
    auto dt1 = E(t,1,deg);

    auto  dptbydt0 = E(pts,0) / d1;
    auto  dptbydt1 = E(pts,1) / d2;

    return vec_t( dptbydt0 -  dt0/dt1 * (dptbydt1 - dptbydt0) );
}

template < class cpts_t, class params_t  >
decltype( pts[0] - pts[1] )
end_vec(const cpts_t & pts,
        const params_t& t, int deg)
{
    return start_vec(rvector(pts),
                     rvector(t), deg);
}

using util::rvector;

template <int dim>
std::vector<pt_t < dim> >&&
hermite_interpolation(const std::vector<double>& t,
                      const std::vector<pt_t<dim> > & pts,
                      vec_t<dim> & v1,
                      vec_t<dim> & v2
                     )
{
    if(t.size() == 0 || t.size()!= pts.size())
        return ;

    using Eigen::TridiagonalMatrix;
    using Eigen::Dynamic;
    TridiagonalMatrix<double,Dynamic> mat(pts.size() + 2);
    rmat_base rm(t, 3);

    mat(0, 0) = rm.der(1, t[0]);
    mat(0, 1) = rm.der(2, t[0]);

    for(size_t j =  1; j <= m; ++j) {
        for(size_t k = 0; k < 3; ++k)
            mat(j, j - 1 + k) = rm.coeff(k, t[j - 1]);
    }

    mat(m + 1,m )    =  rm.der(1, t[m - 1]);
    mat(m + 1,m + 1) =  rm.der(2, t[m - 1]);

    // eigen does not have a banded lu solve...sigh
    banded_lu_decompose(mat);

    std::vector<pt_t<dim> > cpts(pts + 2);

    for(int k = 0; k < dim; ++k) {
        using Eigen::VectorXd;
        VectorXd rhs(m + 2);

        rhs(0) = v1[k];
        for(int l = 1; l < m; ++l)
            rhs(l) = pts[l - 1].p[k];

        rhs(m + 1) = v2[k];
        banded_lu_solve(mat, rhs);

        for(int p = 0; p < m + 2; ++p) {
            cpts[p].p[k] = rhs[p];
        }
    }
    return cpts;
}



}

// {{{ -- interpolation (hermite & global)
std::unique_ptr<Eigen::Map<MatT> >
setUpQMatrix(PointIter pb, PointIter pe)
{
    point_t cg(0);
    cg = centroid(begin(pb), end(pe));
    using boost::make_iterator_range;
    std::array<double,dim*dim> sigmaxy;
    sigmaxy.fill(0.0);

    //should be able to run in parallel
    for( const point_t & p : make_iterator_range(pb,pe))
    {
        for(int i = 0;i < dim;++i)
            for(int j = i ; j < dim;++j)
                sigmaxy[i*dim+j] +=
                    (p[i] - cg[i])*(p[j] - cg[j]);
    }

    for(int i = 0; i < dim;++i) {
        for(int j = i; j < n ; ++j){
            sigmaxy[i][j] /= n;
            sigmaxy[j][i] = sigmaxy[i][j];
        }
    }
    return new Eigen::Map<MatT>(sigmaxy.data());
}

using interpolation_options_t::pametrization_option;
using interpolation_options_t::chord_length;
using interpolation_options_t::centripetal;
using interpolation_options_t::affinely_invariant;
using interpolation_options_t::neilson_foley;
using boost::mpl;

template <class PointIter>
void find_parameters( PointIter pb, PointIter pe,
                      KnotIter tb,
                      mpl::integral_c<parametrization_option,
                      chord_length> )
{
    std::size_t n = std::distance(pb,pe);
    std::vector<double> widths(n);
    widths[0] = 0.0;
    for (int i = 0; i < n - 1; ++i) {
        double dist = 0, l = 0;
        auto v = pr[i + 1] - pr[i];
        l = sqrt(dot(v,v));
        widths[i+1] = l;
    }
    std::partial_sums(widths.begin(),widths.end(),tb);
}

template <class PointIter>
void find_parameters( PointIter pb, PointIter pe,
                      KnotIter tb,
                      mpl::integral_c<parametrization_option,
                      centripetal_length> )
{
    std::size_t n = std::distance(pb,pe);
    std::vector<double> widths(n);
    widths[0] = 0.0;
    for (int i = 0; i < n - 1; ++i) {
        double l = 0;
        auto v = pr[i + 1] - pr[i];
        l = sqrt(sqrt(dot(v,v)));
        widths[i+1] = l;
    }
    std::partial_sums(widths.begin(),widths.end(),tb);
}

template <class PointIter>
void find_parameters(PointIter pb,
                     PointIter pe,
                     KnotIter t,
                     mpl::integral_c<parametrization_option,
                     affinely_invariant>)
{
    std::size_t n = std::distance(pb,pe);
    std::vector<double> widths(n);
    widths[0] = 0.0;
    std::unique_ptr<Map<MatT> > Q( setupQMatrix(pb,pe) );
    for (int i = 0; i < n - 1; ++i) {
        double l = 0;
        auto v = make_eigen_vec ( pr[i + 1] - pr[i] );
        l = sqrt (v.transpose() * (Q.get()) *  v) ;
        widths[i+1] = l;
    }
    std::partial_sums(widths.begin(),widths.end(),t);
}

template <class PointIter>
void find_parameters(PointIter pb, PointIter pe,
                     KnotIter t,
                     mpl::integral_c<parametrization_option,
                     neilson_foley> )
{
    std::size_t n = std::distance(pb,pe);
    std::vector<double> widths(n), thetas(n);
    widths[0] = 0.0;
    thetas[0] = 0.0;
    std::unique_ptr<Map<MatT> > Q( setupQMatrix(pb,pe) );

    for (int i = 0; i < n - 1; ++i) {
        double l = 0;
        auto const & v = make_eigen_vec ( pb[i + 1] - pb[i] );
        decltype(v) u;

        l = sqrt( (v.transpose() * (*Q) *  v)(0) );
        widths[i+1] = l;
        if(i > 0 )
            thetas[i] = angle(u,v);
        u = v;
    }

    for (int i = 1; i < n - 1; ++i)
        width[i] *= (1
                     + 3 *thetas[i] * widths[i - 1]
                     / 2*( widths[i - 1] + widths[i] )
                     + 3 *thetas[i + 1] * widths[i + 1]
                     / 2*( widths[i] + widths[i + 1] ) );

    std::partial_sums(widths.begin(),widths.end(),t);
}

template <class Matrix,class PointIter,
          class KnotIter
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
    int n = std::distance(tb,te);
    BOOST_ASSERT( n > 2);
    const int dim = p[0]->dimension();
    // m(0,1) = E(t,0);
    // m(0,0) = 2*m(0,1);

    for (int j = 1; j < n - 1; ++j) {
        auto aj = E(t,j);
        auto aj1 = E(t,j+1);

        m(j, j - 1) = aj;
        m(j, j)     = 2 * (aj + aj1);
        m(j, j + 1) = aj1;
    }
    // m(n,n-1) = E(t,n-1);
    // m(n,n)   = 2*m(n,n-1);

    for (int j = 1; j < n - 1; ++j)
    {
        double aj = E(t,j);
        double aj1 = E(t,j-1);

        double a = -3 * aj / aj1;
        double b =  3 * (aj / aj1  - aj1 / aj);
        double c =  3 * (aj1 / aj);

        d(j) = a *ev(p[j - 1]) + b*ev(p[j]) +
            c * ev(p[j + 1]);
    }

    // use wherever available, explicit information about derivatives.
    int j = 0;
    for(auto v : vecs)
    {
        if( !tol::eq(sqlen(v), 0) )
        {
            m.row(j).setZero();
            m(j, j) = 1;
            d.row(j) = ev(v);
        }
        ++j;
    }
}

template <class PointIter,class KnotIter,class VecT>
void
eval_tangents_for_pchip( PointIter pb, PointIter pe,
          KnotIter tb, KnotIter te,
          class VecT& explicit_tgts,
          mpl::integral_c<periodic>
    )
{
    static const int dim = point_traits<PointIter>::dimension;
    using Eigen::MatrixXd;
    int n = std::distance(pb,pe);

    if(!tol::eq(pb[0],*pb[n-1],1e-8))
        throw geom_exception(interp_ends_dont_match_for_periodic);

    MatrixXd mat(n);
    Eigen::Matrix<double,Eigen::Dynamic,dim> rhs(n);

    mat(0,0)   = E(tb,n-3);
    mat(0,n-3) = E(tb,n-2);
    mat(0,n-2) = 2*(mat(0,0)  + mat(0,n-3));

    double a0 = (mat(0,n-3)/mat(0,0));

    rhs.row(0) = -3*a0 * ev(E(pb, 0, n-2))
        +3/a0 * ev(E(pb, n-3));

    mat(n-1,n-2) = E(tb,0);
    mat(n-1,1)   = E(tb,n-2);
    mat(n-1,0)   = 2*( mat(n-1,n-2) + mat(n-1,1) );

    a0 = mat(n-1,n-2)/mat(n-1,1);
    rhs.row(n-1) = 3*a0* E(pb,0) - 3/a0 * E(pb,0,n-1);

    setup_mat_for_interior_tgts(mat,pb,pe,tb,te,explicit_tgts,rhs);

    mat.solveInPlace(rhs);

    for(int i = 0; i < n; ++i)
        for(int j = 0; j < dim;++j)
            explicit_tgts[i][k] = rhs[i][k];
}

template <class PointIter,class KnotIter,class VecT>
eval_tangents_for_pchip( PointIter pb, PointIter pe,
          KnotIter tb, KnotIter te,
          VecT& explicit_tgts,
          mpl::integral_c<parabolic_blending> )
{
    int n = std::distance(pb,pe);
    Eigen::TriDiagonalMatrix<double,Dynamic> m(n);
    Eigen::Matrix<double,Eigen::Dynamic,dim> rhs(n);
    m(0,0)     = 1;
    m(n-1,n-1) = 1;
   //see 3.19 expanded hoschek
    rhs.row(0)   = (-3* pb[0] + 4*pb[1] - pb[2])/2;
    rhs.row(n-1) = (pb[n-3] - 4*pb[n-2] + 3*pb[n-1])/2;

    setup_mat_for_interior_tgts(mat, pb, pe, tb, te, explicit_tgts, rhs);
    mat.solveInPlace(rhs);

    for(int i = 0; i < n; ++i)
        for(int j = 0; j < dim;++j)
            explicit_tgts[i][k] = rhs[i][k];
}

template <class PointIter,class KnotIter,class VecT>
eval_tangents_for_pchip( PointIter pb, PointIter pe,
                         KnotIter tb, KnotIter te,
                         VecT& explicit_tgts,
                         mpl::integral_c<vanishing_double_derivatives> )
{
    int n = std::distance(pb,pe);
    Eigen::TriDiagonalMatrix<double,Dynamic> m(n);
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
    a1 = E(params,n-2);
    m(n - 1, n - 2) = a1;
    m(n - 1, n - 1) = 2 * a1;     // there is a mistake here in hoschek

    rhs.row(0) = -3 * pr[0] + 3 * pr[1];
    rhs.row(n - 1) = -3 * pr[n - 2] + 3 * pr[n - 1];


    setup_mat_for_interior_tgts(mat, pb,
                                pe, tb, te,
                                explicit_tgts, rhs);
    mat.solveInPlace(rhs);

    for(int i = 0; i < n; ++i)
        for(int j = 0; j < dim;++j)
            explicit_tgts[i][k] = rhs[i][k];
}

template <class PointIter,class KnotIter,class VecT>
void
eval_tangents_for_pchip( PointIter pb, PointIter pe,
          KnotIter tb, KnotIter te,
          VecT& explicit_tgts,
          mpl::integral_c<not_a_knot> )
{
    Eigen::TriDiagonalMatrix<double,Dynamic> m(n);
    Eigen::Matrix<double,Eigen::Dynamic,dim> rhs(n);
    m(0,0) = (1/a0)* ( 1/a0  +1/a1) ;

    m(0,1) = (1/a0 + 1/a1);// square this
    m(0,1) *= m(0,1);

    auto a0 = E(tb,0);
    auto a1 = E(tb,1);

    auto b0 = 2/(a0 * a0 * a0);
    auto b1 = 3/(a1 * a0 * a0);
    auto b2 = 1/(a1 * a1 * a1);

    rhs.row(0) =
        -(b0 + b1)        * pb[0]
        +(b0  - b2  + b1) * pb[1]
        + b2              * pb[2];

    h0 = E(tb,n-3);
    h1 = E(tb,n-2);

    m(n-1,n-2) = (h0+h1)*(h0+h1)/h1;
    m(n-1,n-1) = (h0 + h0*h0/h1);

    b0 = h0/h1;
    b1 = h1/h0;

    rhs.row(n-1) =
        - b1                  * pb[n-3]
        - (b0 + b1 + 2*b0*b0) * pb[n-2]
        + (2*b0*b0 + 3*b0   ) * pb[n-1];

    setup_mat_for_interior_tgts(mat, pb, pe, tb, te, explicit_tgts, rhs);
    mat.solveInPlace(rhs);

    for(int i = 0; i < n; ++i)
        for(int j = 0; j < dim;++j)
            explicit_tgts[i][k] = rhs[i][k];
}

template <class PointIter, class VectsT,class ParamsT>
std::unique_ptr<bspline< point_traits<PointIter>::point_t > >
pchip_open(PointIter pb, PointIter pe, const VecsT& tgts,
           const ParamsT& params)
{
    static const int dim = point_traits<PointIter>::dimension;
    std::vector<pt_t<dim> > cpts;
    std::vector<double> knots;
    cpts.resize( 2* n );
    knots.resize( 2* n + 4);

    cpts[0] = pb[0];

    cpts[1] = pb[0] + 1/3 * tgts[0] * E(params,0);

    knots[0] = knots[1] = knots[2] = knots[3] = params[0];

    for( int i = 2; i < n ; ++i) {
        cpts[2*i-2] = pb[i-1] - 1/3 * tgts[i-1] *  E(params,(i-2)) ;

        cpts[2*i-1] = pb[i-1] + 1/3 * tgts[i-1] *  E(params,(i-1)) ;

        knots[2*i] = knots[2*i+1] = params[i];
    }

    cpts[2*n-2] = pb[n-1] - 1/3 * tgts[n-1] * E(params,n-2);
    cpts[2*n-1] = pb[n-1];
    knots[2*n+3] = knots[2*n+2] = knots[2*n+1] = knots[2*n] = params[n-1];
    return std::unique_ptr<bspline<Point> >(
        new bspline<Point>(cpts, knots, 3));
}

template <class PointIter, class VectsT,class ParamsT>
std::unique_ptr<bspline< point_traits<PointIter>::point_t > >
pchip_closed(PointIter pb, PointIter pe,
             const VecsT& tgts,
             const ParamsT& params)
{
    cpts.resize( 2* n );
    knots.resize( 2* n );
//pg 180,188 hoschek
    cpts[0] = pr[0] - 1/3 * row(d,0) *  E(params,n-2) ;

    cpts[1] = pr[0] + 1/3 * row(d,0) *  E(params, 0) ;

    knots[0] = knots[1] = params[0];

    for( int i = 2; i < n ; ++i ) {
        cpts[ (2*i-2) ] = pr[i-1]- 1/3 *row(d,i-1) * E(params,(i-2));

        cpts[ (2*i-1) ] = pr[i-1] + 1/3 * row(d,i-1) * E(params,(i-1)) ;

        knots[2*i-2] = knots[2*i-1] = params[i];
    }

    cpts[2*n-2] = pr[n-1] - 1/3 * row(d,n-1) * E(params,n-2) ;

    cpts[2*n-1] = pr[n-1] + 1/3 * row(d,n-1) * E(params,0);

    knots[2*n-2] = knots[2*n-1] = params[n-1];

    BOOST_ASSERT( pr[n-1] == pr[0]);

    return std::unique_ptr<bspline< Point> > (
        new bspline<Point>( cpts, knots , 3));
}

template <class PointIter, class VecT>
std::unique_ptr< bspline< dim  > >
piecewise_cubic_hermite_interp
(
    PointIter pb,
    PointIter pe,
    interpolation_options_t opts,
    const VecT & vecs
)
{
    using  point_t = point_traits<PointIter>::point_t;
    typedef bspline<point_t>   bspline_t;
    static const int dim =  point_traits<point_t>::dim;

    // generate parameters using chord length approximation.
    const int m = std::distance(pb,pe);
    BOOST_ASSERT( opts.end_conditions !=
                  interpolation_options_t::periodic
                  || pb[n-1] == pb[0] );
    BOOST_ASSERT( n + vecs.size() >=4);

    std::vector < double > params(m);

    switch_case_4(opts.parametrization,
                  centripetal,
                  find_parameters,
                  pb, pe, params.begin());


    typedef Eigen::Matrix<double,Dynamic,dim> MatrixXdim;
    typedef Eigen::Matrix<double,Dynamic,Dynamic> MatrixXX;

    MatrixXdim tgts(m, dim);
    MatrixXX   mat(m, m); // to determine tangents

// end conditions depend on options passed..
    bool periodic = opts.end_conditions == periodic;

    switch_case5( opts.end_condition, periodic,
                  eval_tangents_for_pchip,err,pb,pe,tb,tw,tgts);

    std::vector<point_t> cpts;
    std::vector<double> knots;
    smooth_tgts(params,tgts);
    if( !periodic) {
        return pchip_open(tb,te,params,tgts);
    }else {
        return pchip_closed(tb,te,params,tgts);
    }
}

// fill knots
template < class KnotIter, class BoolRange , class OutputKnotIter>
void fill_knots (int d, KnotIter param_f, KnotIter param_l,
                 const  BoolRange &hasTangents, OutputKnotIter t)
{
    int k = 0;

    t = std::fill_n(t, d + 1, *param_f);

    double last_avg = 0;
    int j =0;
    for (; param_f != param_l; ++param_f, ++j)
    {
        double avg = std::accumulate(param_f, param_f + d, 0) / d;

        if (hasTangents[j])
        {
            *t++ = (avg + last_avg) / 2;
        }
        *t++ = avg;
        last_avg = avg;
    }
    std::fill_n(t, d + 1, *(param_l - 1));
}
/*
// interpolation using point tange and param range
template < class PointRange, class ParamRange, class VectorRange >
std::unique_ptr < bspline< boost::range_value<PointRange>::type::dim > >
global_curve_interpolate(const PointRange & point_range,
                         const VectorRange & tangent_range,
                         const ParamRange & param_range,
                         int degree, bool periodic)
{

    BOOST_ASSERT( degree >=2);
    BOOST_ASSERT(size(param_range) == size(point_range));
    BOOST_ASSERT(size(tangent_range) >= size(point_range));

    auto tangent_f = begin(tangent_range);
    auto tangent_l = end(tangent_range);

    const int d = degree;
    int n = boost::size(point_range);

    typedef vec_t point_t;
    typedef bspline<point_t > bspline_t;
    typedef periodic_bspline<point_t > periodic_bspline_t;

    //   typedef bspline<typename boost::range_value<PointRange>::type > bspline_t;

    if( n <2) {
//TODO: throw an exception
        return std::unique_ptr<bspline_t>();
    }

    std::unique_pointer<bool[]> hasTangents(new bool[n]);

    const int dim = point_traits<vec_t>::dim;

    int num_tgts = 0;
    int j =0;
    auto tgt_f = tangent_f;
    for(;tgt_f != tangent_l;++tgt_f ) {
        vec_t p = *tgt_f;
        bool tangentProvided  =  ( p != point_traits<vec_t>::zero_pt() );
        if( tangentProvided ) {++num_tgts;}
        hasTangents[j++] = tangentProvided;
    }

    int num_cpts = n + num_tgts;
    std::unique_array < double > t( new double[num_cpts + d + 1]);
    int num_knots = num_cpts + d +1;

    // todo: for periodic splines.
    auto param_f = begin(param_range);
    auto param_l = end(param_range);

    // fill up knots
    fill_knots(d,param_f, param_l, hasTangents, t.get());

    //{{{
    // set up linear system by fill constructinng the banded matrix m
    // pg 116 of spline methods (tom lyche)
    banded_matrix < double > m(num_cpts, num_cpts, d + 1, d + 1);

    int k = 0;
    tgt_f = tangent_f;
    vector<double> coeffs[2];
    for (int i = 0; i < num_knots; ++i , ++tgt_f)
    {
        basis_derivatives_at(param_range[i],
                             t.get(), t.get() + num_knots,
                             d, 0, 1,coeffs);

        std::copy_n( coeffs[0].begin(), d+1, m.find1(0,k++,i));

        if( hasTangents[i] )
        {
            // a tangent is provided, so we have et another linear eqn.
            std::copy_n( coeffs[1].begin(), d+1, m.find1(0,k++,i));
        }
    }
    //}}}

    auto point_f = begin(point_range);
    auto point_l = end( point_range);
    tgt_f = tangent_f;


    //{{{  set up b of the linear system (mx =b) consisting of
    //    points and tangent (if present)
    vector< vec_t > b(num_cpts);
    j =0;
    for(int i=0 ; i < n;++tgt_f,++i)
    {
        b[j++] = point_f[i];
        if( hasTangents[i])
        {
             for( int k =0; k < dim; ++k) {
                 b[j][k] = (*tgt_f)[k];
             }
             ct_copy<dim>::do_copy(  *tgt_f, b[j] );
             ++j;
        }
    }
    //}}}
    lu_solve(m,b); // banded lu solve

    if( periodic ) {
        return std::unique_ptr < bspline_t >(
            new periobdic_bspline_t (b,   boost::make_iterator_range(t.get(), t.get() + num_knots), d ) );
    }
    return std::unique_ptr < bspline_t >(
        new bspline_t (b,
                       boost::make_iterator_range(t.get(), t.get() + num_knots),
                       d ) );
}
*/
// }}}
// }}}

#include "util.hpp"
#include "Eigen/Core"
#include <Eigen/Dense>
#include "point.hpp"
#include "bspline.hpp"
#include "periodic_bspline_cons.hpp"
#include <type_traits>
#include <numeric>
#include "geom_exception.hpp"
#include "tol.hpp"

namespace geom {
template < class cpts_t, class param_t  >
auto
start_vec(const cpts_t & pts,
          const param_t& t, int deg)->decltype( pts[0] - pts[1] )
{
    typdef decltype( pts[0] - pts[1] ) vec_t;
    auto dt0 = E(t,0,deg);
    auto dt1 = E(t,1,deg);

    auto  dptbydt0 = E(pts,0) / d1;
    auto  dptbydt1 = E(pts,1) / d2;

    return vec_t( dptbydt0 -  dt0/dt1 * (dptbydt1 - dptbydt0) );
}

template < class cpts_t, class params_t  >
auto
end_vec(const cpts_t & pts,
        const params_t& t, int deg)->decltype( pts[0] - pts[1] )
{
    return start_vec(rvector(pts),
                     rvector(t), deg);
}

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

    point_iter_traits<pt_t<dim>*>::PointContT cpts(pts + 2);

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

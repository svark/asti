// -*- mode:c++ -*-
#include "box_compute.hpp"
#include <math.h>
#include <utility>
#include <Eigen/Core>
#include <Eigen/SVD>
#include "tol.hpp"
#include "point.hpp"
#include "reparametrize.hpp"
#include "implicit.hpp"
#include "any_curve.hpp"
#include "split_into_bezier_patches.hpp"
// find a basis q_k of M implicit polynomials
// expand q_k(p(t))  in terms of \alpha_k(t)
// coefficients of \alpha_k(t) will form matrix D
// q(p(t))= (Db)^T [\alpha_k(t)] find b of norm 1
// corresponding to min singular value

namespace geom {

//{{{  (@* "coordinate transformed to lie within unit box")
template <class Point>
struct homogc
{
    template <class Curve>
    homogc(Curve c)
    :s(make_any_curve(std::move(c)))
    {
    }
	homogc(homogc&& o)
		:s(std::forward<any_curve < Point > >(o.s)){}

	homogc(const homogc& o)
		:s(o.s){}


    Point eval(double param) const
    {
        return s.eval(param);
    }

    Point operator()(const Point& p) const
    {
        return p;
    }
private:
    any_curve <Point>  s;
};

double
implicitCurveFormBase::eval(const point2d_t& pt) const
{
	return eval(point3d_t(pt,1.0));
}

double
implicitCurveFormBase::eval3d(const point3d_t& p3d) const
{
	long c0  = 1;
    int j = 0;
    double qsum = 0;


    for(int ks = 0; ks <= qdeg; ++ks)
    {
        int k3 = qdeg - ks;
        if(ks) {
            c0 *= (qdeg - ks + 1);
            c0 /= ks;
        }

        long c1 = 1;
        for(int k1 = 0; k1 <= ks; ++k1,++j)
        {
            int k2 = ks - k1;
            if(k1) {
                c1 *= (ks - k1 + 1);
                c1 /= k1;
            }
            auto c = c0 * c1;  //c == (k1+k2+k3)!/(k1!k2!k3!)

            qsum += coefficients[j] * c * std::pow(p3d[0], k1)
                *std::pow(p3d[1], k2)
                *std::pow(p3d[2], k3);
        }
    }
    return qsum;
}
//}}}
typedef homogc<point3d_t> homogc3d_t;
typedef homogc<point2d_t> homogc2d_t;
//{{{  (@* "2d and 3d implicit forms")
struct implicitCurveForm3d : public implicitCurveFormBase
{
    implicitCurveForm3d(std::vector<double> coefficients_,
                        int qdeg_,
                        homogc3d_t hg_)
        :implicitCurveFormBase
		(std::move(coefficients_), std::move(qdeg_)), hg(std::move(hg_))
    {
    }
    double eval(const point3d_t & p) const
    {
        return eval3d(hg(p));
    }
private:
  homogc3d_t hg;
};

//}}}

//{{{  (@* "compute eigen vector correcponding to min eigen value")
void
computeSingularVec(const Eigen::MatrixXd & mat, std::vector<double>& vs)
{
	Eigen::JacobiSVD<Eigen::MatrixXd> svd;
	svd.compute(mat,Eigen::ComputeFullV);


    auto& evs =  svd.matrixV();
    auto& v =  evs.col(evs.cols()-1) ;
	int sz = v.size();

    for(int i = 0; i < sz; ++i)
    {
		vs.push_back(v[i]);
    }
}
//}}}

// (@file :file-name "./media/implicit.pdf" :to "./media/implicit.pdf" :display "implicitize")
//{{{ (@* "implicitize a homog curve")
std::unique_ptr<implicitCurveFormBase>
implicitize( const homogc < point3d_t >& hg, int qdeg, int sdeg)
{
	long c0 = 1;
    long L = qdeg*sdeg + 1;       // number of basis functions \alpha_k(t)
    long M = (qdeg + 2) * (qdeg + 1) / 2; // number of monomial basis functions
    Eigen::MatrixXd mat(L,M);
    int  j = 0;

	std::vector<double> ks(2 * L);
	std::fill_n(ks.begin(), L, 0);
	std::fill_n(ks.begin()+L, L, 1.0);

    int basisdeg = L - 1;
    rmat_base_vd rm(ks, basisdeg);
    for(int k12 = 0; k12 <= qdeg; ++k12)
    {
        int k3 = qdeg - k12;
        if(k12) {
            c0 *= (qdeg - k12 + 1);
            c0 /= k12;
        }

        long c1 = 1;

        for(int k1 = 0; k1 <= k12; ++k1,++j)
        {
            int k2 = k12 - k1;
            if(k1) {
                c1 *= (k12 - k1 + 1);
                c1 /= k1;
            }
            auto c = c0 * c1;
            //c = (k1+k2+k3)!/(k1!k2!k3!)
            //  =  (k1+k2+k3)!/(k1+k2)!k3! \times (k1+k2)!/k1!k2!
            Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic, Eigen::RowMajor> qmat(L, L);
            qmat.setZero();
            //  determine c_k such that q_k(\gamma(t)) =  sum{ c_k
            //  \alpha_k(t) }
            //  {\alpha_k(t)} are bernstein basis functions of degree L - 1
            Eigen::VectorXd rhs(L);
            for(int i = 0;i < L; ++i) {
                double  u =   double(i) / (basisdeg);
                std::vector<double>  basis;

                std::tie(basis, std::ignore) = rm.get_basis(u);
                for(int k = 0;k < L; ++k)  {
                    qmat(i, k) = basis[k];
                }
                auto p3d = hg.eval(u);
                rhs(i) =
                    c * std::pow(p3d[0],k1)
                      * std::pow(p3d[1],k2)
                      * std::pow(p3d[2],k3);
            }
            mat.col(j) = qmat.lu().solve(rhs);
#ifndef NDEBUG
			std::vector<double> cpts(L);
			for(int l = 0;l  < L;++l)
			{
				cpts[l] = mat(l,j);
			}
			rmat<double> rm_chk(cpts, ks,basisdeg);
			assert( tol::eq(rm_chk.eval(0) , rhs(0)) );
			assert(	tol::eq(rm_chk.eval(1.0/basisdeg), rhs(1)) );
			assert(	tol::eq(rm_chk.eval(2.0/basisdeg) , rhs(2)) );
			assert(	tol::eq(rm_chk.eval(1), rhs(L-1)) );
#endif
        }
    }

    std::vector<double> vs;
    vs.reserve(M);
    computeSingularVec(mat, vs);

    return
        std::unique_ptr < implicitCurveFormBase >
        (new implicitCurveForm3d(vs, qdeg, hg));
}

//}}}
//{{{  (@* "implicitize 2d rational bspline")
std::vector<std::unique_ptr < implicitCurveFormBase> >
implicitize(const rational_bspline<point2d_t>& spl, int qdeg)
{
    const int sdeg     = spl.degree();

	auto const & patches = ops::split_into_bezier_patches(spl);
	std::vector<  std::unique_ptr<implicitCurveFormBase> > imps;
	for( auto const &s:  patches ) {
	   homogc < point3d_t > hg(ops::reparametrize(s.spline(),0,1));
	   imps.emplace_back(implicitize(hg,qdeg,sdeg) );
	}
	return imps;
}

std::unique_ptr<implicitCurveFormBase>
implicitize( const std::function<point3d_t(double)>& f, int qdeg,int sdeg)
{
	homogc < point3d_t > hg(f);
	return implicitize(hg, qdeg, sdeg);
}
//}}}

//{{{  (@* "implicitize 2d bspline")
std::vector<std::unique_ptr < implicitCurveFormBase> >
implicitize(const bspline<point2d_t>& spl, int qdeg)
{
   auto &cpts = spl.control_points();
   std::vector<double> weights(cpts.size(), 1);
   auto rbs (make_rbspline(cpts, weights, spl.knots(), spl.degree()));
   return implicitize(rbs, qdeg);

}
//}}}

}

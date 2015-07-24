// -*- mode:c++ -*-
#ifndef ASTI_IMPLICIT_HPP
#define ASTI_IMPLICIT_HPP
#include "box_compute.hpp"
#include <math.h>
#include <utility>
#include <Eigen/Core>
#include <Eigen/EigenValues>
#include "tol.hpp"
#include "point.hpp"
#include "reparametrize.hpp"
#include "implicit.hpp"

// find a basis q_k of M implicit polynomials
// expand q_k(p(t))  in terms of \alpha_k(t)
// coefficients of \alpha_k(t) will form matrix D
// q(p(t))= (Db)^T [\alpha_k(t)] find b of norm 1
// corresponding to min eigen value

namespace geom {

//{{{  (@* "coordinate transformed to lie within unit box")
template <class Point>
struct homogc
{
    typedef RAWTYPE(make_vec(Point())) VectorT;
    homogc(const bspline <Point> &s_)
        :s(s_),center(new VectorT(0.0))
    {
        auto  b  = geom::ops::compute_box(s);
        halfdiag = len(make_vec(b.size(X) / 2, b.size(Y) / 2));
        *center  =  make_vec( lerp(0.5, b.lo, b.hi) );
    }
	homogc(homogc&& o)
		:s(std::forward<bspline<Point> >(o.s)),halfdiag(o.halfdiag),
		 center(std::move(o.center)){}

	homogc(const homogc& o)
		:s(o.s),halfdiag(o.halfdiag),
		 center(new VectorT(*o.center)){}

	homogc(bspline <Point> &&s_)
        :s(std::forward<bspline<Point>>(s_)),center(new VectorT(0.0))
    {
        auto  b  = geom::ops::compute_box(s);
        halfdiag = len(make_vec(b.size(X) / 2, b.size(Y) / 2));
        *center   =  make_vec(lerp(0.5, b.lo, b.hi));
	}

    Point eval(double param) const
    {
        return s.eval(param);
    }

    Point operator()(const Point& p) const
    {
        return p;
    }
private:
    bspline <Point>  s;
    double           halfdiag;

	std::unique_ptr<VectorT>    center; //made it a pointer to prevent alignment issues
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
computeEigenVec(const Eigen::MatrixXd & mat, std::vector<double>& vs)
{
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(mat, Eigen::ComputeThinU|Eigen::ComputeThinV);
	//Eigen::VectorXd evals = svd.singularValues();
	
    int       minIndex     = svd.nonzeroSingularValues()-1;

    assert(minIndex!=- 1);
    // std::vector<double> vs;
    auto& evs =  svd.matrixV();
    auto& v =  evs.col(minIndex) ;
	int rows = v.size();
    for(int i = 0; i < rows; ++i)
    {
		vs.push_back(v[i]);
    }
}
//}}}

// (@file :file-name "./media/implicit.pdf" :to "./media/implicit.pdf" :display "implicitize")

//{{{  (@* "implicitize 2d rational bspline")
std::unique_ptr < implicitCurveFormBase >
implicitize(const rational_bspline<point2d_t>& spl, int qdeg)
{
    const int sdeg     = spl.degree();
    homogc < point3d_t > hg( ops::reparametrize(spl.spline(),0,1));

    long c0 = 1;
    long L = qdeg*sdeg + 1;       // number of basis functions \alpha_k(t)
    long M = (qdeg + 2) * (qdeg + 1) / 2; // number of monomial basis functions
    Eigen::MatrixXd mat(L,M);
    int  j = 0;

	std::vector<double> ks(2 * L, 0);
    for(unsigned long i = L;i < ks.size(); ++i)
        ks[i] = 1.0;
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
            Eigen::MatrixXd qmat(L, L);
            qmat.setZero();
           
            Eigen::VectorXd rhs(L);
            for(int i = 0;i < L; ++i) {
                double  u =  double(i) / (basisdeg);
                std::vector<double>  basis;
                int nu;
                std::tie(basis, nu) = rm.get_basis(u);
                for(int k = 0;k < L; ++k)  {
                    qmat(i, nu - basisdeg + k) = basis[k];
                }
                auto p3d = hg.eval(u);
                rhs(i) =
                    c* std::pow(p3d[0],k1)
                    * std::pow(p3d[1],k2) * std::pow(p3d[2],k3);
            }
            mat.col(j) = qmat.lu().solve(rhs);
        }
    }

    std::vector<double> vs;
    computeEigenVec(mat, vs);

    return
        std::unique_ptr < implicitCurveFormBase >
        (new implicitCurveForm3d(vs, qdeg, hg));
}
//}}}

//{{{  (@* "implicitize 2d bspline")
std::unique_ptr < implicitCurveFormBase > 
implicitize(const bspline<point2d_t>& spl, int qdeg)
{
   auto &cpts = spl.control_points();
   std::vector<double> weights(cpts.size(), 1);
   auto rbs (make_rbspline(cpts, weights, spl.knots(), spl.degree()));
   return implicitize(rbs, qdeg);
   
}
//}}}

}
#endif //ASTI_IMPLICIT_HPP

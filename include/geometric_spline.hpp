#ifndef ASTI_GEOMETRIC_SPLINE
#define ASTI_GEOMETRIC_SPLINE
#include "split_into_bezier_patches.hpp"
namespace geom
{

//
// Suppose two parametric curves $X(u)$ for $u \in [u_0,u_1]$ and $Y(t)$
// for $t \in [t_0,t_1]$, lying  in %\R^d$ join continuously at a
// regular point $P = X(u)1) = Y(t_0)$ (i.e derivatives are not zero)
// then we say that they join with $C^1$ continuity provided that at $P$
// the derivatives are the same . On the other hand, the weaker more
// geometric condition that the tangents have the same direction (but
// possibly different lengths) is already satisifed when the derivative
// vectors point in the same direction. Thus tangent continuity is
// written as
//\[
//  Y'(t_0^+) = \omega_11 X'(u_1^{-})
//\]
// where $omega_11 > 0$ in view of the regularity assumption.
// this is called GC^1 continuity or FC^1 (frenet frame)
// continuity. (Note FC^r may not be same as GC^r for r > 2)
// from the equations for computing the curvature and torsion of the
// curve, we sae that provided the FC^1 condition holds, curvature
// continuity is equivalent to
//\[
//      Y^{''}(t_0^+) = \omega_{11}^2 X^{''}(u_1^{-1}) + \omega_{12} X'(u_1^{-1})
//\]

// for FC^3 continuity (torsion continuity) we need
//
//\[
// Y^{'''}(t_0^{+}) = \omega_{11}^3 X^{'''}(u_1^{-1}) + \omega_{12}
//  X^{'''}(u_1^{-1} + \omega_{23} X'(u_1^{-1})
//\]

//  GC^3 continiuity is a spcial form of torsion continuous splines for which we require
// $\omega{12} = 3 \omega_{23} * \omega_{11}
//

// hoshek pg 222
template <class SplineCurve>
std::pair<SplineCurve,SplineCurve>
make_tangent_continuous(const SplineCurve& s1,
                        const SplineCurve& s2,
                        double omega11/*design param*/)
{
    int deg1 = s1.degree();
    int deg2 = s2.degree();

    typedef typename SplineCurve::cpts_t cpts_t;
    SplineCurve spl1(last_bezier_patch(s1));
    cpts_t b    = spl1.control_points();
    auto   tr2  = s2.knots().begin();
    auto   tr1  = spl1.knots().begin();
    double q    = tr2[deg2 + 1] - tr2[deg2];
    q          /= (tr1[deg1 +1] - tr1[deg1]);
    double N1   = double(deg1)/deg2;
    auto tweak2pts = [&b,&q,&N,&omega11]
        (typename cpts_t::iterator btilde)
        {
            btilde[0] = b.back();
            auto br   = b.rbegin();
            btilde[1] = b.back()
            + (q*N1*omega11)*(br[0] - br[1]);
        };
    return std::make_pair(s1,transform_at_left(s2, tweak2pts));
}

template <class SplineCurve>
std::pair<SplineCurve,SplineCurve>
make_curvature_continuous(const SplineCurve& s1,
                          const SplineCurve& s2,
                          double omega11/*design param*/,
                          double mu /*design param*/)
{
    int deg1 = s1.degree();
    int deg2 = s2.degree();

    typedef typename SplineCurve::cpts_t cpts_t;
    SplineCurve spl1(last_bezier_patch(s1));
    cpts_t      b      = spl1.control_points();
    auto        tr2    = s2.knots().begin();
    auto        tr1    = spl1.knots().begin();
    double      q      = tr2[deg2 + 1] - tr2[deg2];
    q                 /= (tr1[deg1 +1] - tr1[deg1]);
    double      N1     = double(deg1)/deg2;
    double      N2     = deg1*(deg1-1.0)/(deg2*(deg2-1.0));
    auto        tweak3pts  = [&b,&q,&N,&omega11]
        (typename cpts_t::iterator btilde)
        {
            btilde[0] = b.back();
            auto br = b.rbegin();
            btilde[1] = b.back()
            + (q*N1*omega11)*(br[0] - br[1]);
            btilde[2]  = b[1];
            btilde[2] -= q*q*N2*omega11*omega11 * (br[1] - br[2]);
            btilde[2] += mu * (br[0] - br[1]);
        };
    return std::make_pair(s1,transform_at_left(s2, tweak3pts));
}
}
}
#endif // ASTI_GEOMETRIC_SPLINE

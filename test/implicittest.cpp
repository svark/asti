#include <vector>
#include "implicit.hpp"
#include "catch.hpp"
#include "testutils.hpp"
#include "bspline.hpp"
#include "point.hpp"
#include "circle.hpp"
#define _USE_MATH_DEFINES
#include <math.h>
TEST_CASE("implicit2d",  "[bspline][approximate][implicitization]") {
  /*
  cpts2d_t pts(5);
  pts[0] = make_pt(0.0,0.0);
  pts[1] = make_pt(1.0,1.0);
  pts[2] = make_pt(1.5,0.3);
  pts[3] = make_pt(1.8,0.1);
  pts[4] = make_pt(2.0,0.0);
  std::vector<double> ks(3+3+2);
  ks[0] = ks[1] = ks[2] = 0.0;
  ks[3] = 0.5; ks[4] = 0.8;
  ks[5] = ks[6] = ks[7] = 1.0;
  geom::bspline<geom::point2d_t> bs(std::move(pts), std::move(ks),2);
  std::unique_ptr < geom::implicitCurveFormBase> icf = geom::implicitize(bs, 6);*/

  
  using geom::make_pt;
  const double sqrt2 = M_SQRT1_2;
  //auto circ = geom::make_circle(make_pt(1,0), make_pt(sqrt2,sqrt2),make_pt(0,1));
  //auto bsc = to_rational(circ);
  std::function<point3d_t(double )> f = [](double t)  { return make_pt(2*t, 1-t*t,1+t*t); };

  auto icf = geom::implicitize(f, 2, 2);
  REQUIRE( tol::eq(icf->eval(f(1/2.0)), 0, 1e-2));
  REQUIRE( tol::eq(icf->eval(f(1/3.0)), 0, 1e-2));
  REQUIRE( tol::eq(icf->eval(f(1/6.0)), 0, 1.1e-2));
  // REQUIRE(icf.eval(bs.eval(0.5)), Approx(0.0));
  // REQUIRE(icf.eval(bs.eval(0)),  Approx(0.0));
  // REQUIRE(icf.eval(bs.eval(1)), Approx(0.0) );
}

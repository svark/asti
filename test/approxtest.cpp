#include "point.hpp"
#include <vector>
#include "catch.hpp"
#include "testutils.hpp"
#include "spline_approx.hpp"
#include <functional>

TEST_CASE("approxtestbasic",  "[bspline][approx]") {

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

  const double ts[] = {-1,-1,-1,-1,0.5,1,1,1,1};
  auto cube = [](double x) { return x*x*x -3*x+1;};
  auto spl0 = geom::ops::cubic_approx1d(
      std::function<double(double)>(cube),
      std::vector<double>(ts, ts + sizeof(ts) / sizeof(ts[0])));
  INFO( "cube approx at -1,0.1,0,1:"
            << spl0.eval(-1) << ","
            << spl0.eval(0.1) << ',' << spl0.eval(0) << ','
        << spl0.eval(1) << "\n" );
  REQUIRE(spl0.eval( - 1) == Approx(3));
  REQUIRE(spl0.eval( 0.1) == Approx(0.701));
  REQUIRE(spl0.eval( 0)   == Approx(1));
  REQUIRE(spl0.eval( 1)   == Approx( - 1));
  auto p = make_pt(0.68,0.586);
  std::vector<double> newknots(bs.knots().size() + 2);
  double arr[] = {0.225,0.75};
  std::merge(bs.knots().cbegin(), bs.knots().cend(), arr, arr + 2, newknots.begin());
  newknots.insert(newknots.begin(), 0.0);
  newknots.insert(newknots.end(), 1.0);
  auto v_vdash = [&p,&bs]
    (double u) -> double
      {
    auto v =   (bs.eval(u) - p);
    auto vdash =  bs.eval_derivative(1,u);
    return dot(v, vdash);
  };

  auto spl2 = geom::ops::cubic_approx1d(
      std::function<double(double)>(v_vdash), newknots);
  auto pt_ = spl2.eval(-1);
    INFO(  "\n dist.dist'  at 0.1,0.2,1:" << v_vdash(0.25)
         << ',' << v_vdash(0.2) << ',' << v_vdash(.5) << "\n");
  INFO(  "\n dist.dist' approx at 0.1,0.2,1:"
             << spl2.eval(0.25) << ',' << spl2.eval(0.2) << ','
         << spl2.eval(.5) << "\n");

  REQUIRE( fabs(v_vdash(0.2) - spl2.eval(0.2)) < 0.05);
  REQUIRE( fabs(v_vdash(0.2) - spl2.eval(0.2)) < 0.05);
  REQUIRE( fabs(v_vdash(0.5) - spl2.eval(0.5)) < 0.1);

  auto fp =  geom::ops::foot_param(bs,p);
  INFO( "soln NEAR param .2::" << fp);
  REQUIRE(fp ==  Approx(0.199843));
}

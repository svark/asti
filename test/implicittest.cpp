#include <vector>
#include "implicit.hpp"
#include "catch.hpp"
#include "testutils.hpp"
#include "bspline.hpp"
#include "point.hpp"

TEST_CASE("implicit2d",  "[bspline][approximate][implicitization]") {
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
  std::unique_ptr < geom::implicitCurveFormBase> icf = geom::implicitize(bs, 4);
  REQUIRE(icf->eval(bs.eval(0.2)) == Approx(0.0));
  // REQUIRE(icf.eval(bs.eval(0.5)), Approx(0.0));
  // REQUIRE(icf.eval(bs.eval(0)),  Approx(0.0));
  // REQUIRE(icf.eval(bs.eval(1)), Approx(0.0) );
}

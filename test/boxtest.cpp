//-*- mode:c++ -*-
#include "point.hpp"
#include <vector>
#include "catch.hpp"
#include "testutils.hpp"
#include "box_compute.hpp"

#include "bspline.hpp"
#include "conic.hpp"
#include "tessellate.hpp"

TEST_CASE("boxtest", "[bspline][box]") {
  using geom::make_pt;
  cpts2d_t pts(5);
  pts[0] = make_pt(0.0, 0.0);
  pts[1] = make_pt(1.0, 1.0);
  pts[2] = make_pt(1.5, 0.3);
  pts[3] = make_pt(1.8, 0.1);
  pts[4] = make_pt(2.0, 0.0);
  std::vector<double> ks(3 + 3 + 2);
  ks[0] = ks[1] = ks[2] = 0.0;
  ks[3] = 0.5;
  ks[4] = 0.5;
  ks[5] = ks[6] = ks[7] = 1.0;
  geom::bspline<geom::point2d_t> bs(std::move(pts), std::move(ks), 2);
  using namespace geom::ops;
  auto const &bx = compute_box(bs);
  REQUIRE(bx.lo == make_pt(0, 0));
  REQUIRE(bx.hi == make_pt(2, 1));

  auto const &bxtight = compute_box_tight(bs);
  REQUIRE(bxtight.lo == make_pt(0, 0));
  REQUIRE(bxtight.hi == make_pt(2, 0.588235294));
  auto pl = geom::tessellate(bs, 0.001);
  geom::box<geom::point2d_t>  b;
  for(auto cp : pl.control_points())
  {
	  REQUIRE(geom::contains_point( bxtight, cp, 0.001) );
  }

  double ps[][2] = {{0.,0},{0.0833333,0.5},{0.333333,1.25},{0.666667,0.75},{0.916667,0.5},{1.,0}};
  pts.resize(6);
  for(int i = 0 ; i < 6;++i)
      pts[i] = make_pt(ps[i][0],ps[i][1]);
  ks.resize(6+3+1);
  double karr[] = {0,0,0,0,1.0/4,3.0/4,1,1,1,1};
  for(int i =0; i<sizeof(karr)/sizeof(double);++i)
      ks[i] = karr[i];

  geom::bspline<geom::point2d_t> bs2(std::move(pts),
                                     std::move(ks), 3);
  auto const &bxtight2 = compute_box_tight(bs2);
  REQUIRE(bxtight2.lo == make_pt(0, 0));
  REQUIRE(bxtight2.hi == make_pt(1,0.98));
}

TEST_CASE("rational_boxtest", "[rational][bspline][box]") {
  using geom::make_pt;
  geom::point2d_t pts[3];
  pts[0] = make_pt(0.0, 0.0);
  pts[1] = make_pt(1.0, 1.0);
  pts[2] = make_pt(2.0, 0.5);

  auto ca = geom::make_circular_arc(pts);
  INFO("e:" << ca.eval(0) << " : " << ca.eval(0.5) << " : " << ca.eval(1.0) );
  auto const & rbs  = geom::make_rbspline_from_conic(ca);
  INFO("e:" << rbs.eval(0) << " : " << rbs.eval(0.5) << " : " << rbs.eval(1.0) );
  using namespace geom::ops;
  auto const &bx = compute_box_tight(rbs);
  REQUIRE(bx.lo == make_pt(0, 0));
  REQUIRE(bx.hi == make_pt(2, 1.04020148));

}

//-*- mode:c++ -*-
#include "point.hpp"
#include <vector>
#include "catch.hpp"
#include "testutils.hpp"
#include <functional>
#include "bspline.hpp"
#include "find_param.hpp"
TEST_CASE("find_param",  "[bspline][point inversion]") {

    cpts2d_t pts(5);
    pts[0] = make_pt(0.0,0.0);
    pts[1] = make_pt(1.0,1.0);
    pts[2] = make_pt(1.5,0.3);
    pts[3] = make_pt(1.8,0.1);
    pts[4] = make_pt(2.0,0.0);
    std::vector<double> ks(3+3+2);
    ks[0] = ks[1] = ks[2] = 0.0;
    ks[3] = 0.5; ks[4] = 0.8;
    ks[5] = ks[6] = ks[7] = 1.2;
    geom::bspline<geom::point2d_t> bs(std::move(pts), std::move(ks),2);
    bool ok;
    double u;
    std::tie(u, ok) = geom::find_param(bs.eval(0.1), bs);
    REQUIRE(ok);
    REQUIRE(u == Approx(0.1));
    std::tie(u, ok) = geom::find_param(bs.eval(0.1) + make_vec(0, 0.1), bs);
    REQUIRE(!ok);
}

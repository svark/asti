#include "point.hpp"
#include <vector>
#include "catch.hpp"
#include "bspline.hpp"
#include "testutils.hpp"
#include "periodic_bspline_cons.hpp"
#include "raise_degree.hpp"

TEST_CASE("raise_degree", "[periodic_bspline][degree]")
{
    cpts2d_t pts(4);
    pts[0] = make_pt(0.0,0.0);
    pts[1] = make_pt(0.4,0.3);
    pts[2] = make_pt(0.2,0.8);
    pts[3] = make_pt(-0.2, 0.4);
    //pts[4] = make_pt(0.0,0.0);
    std::vector<double> ks(5);
    ks[0] = 0.0; ks[1] = 0.3; ks[2] = 0.6;
    ks[3] = 0.8; ks[4] = 1.0;
    geom::periodic_bspline<geom::point2d_t> bs( geom::make_periodic_bspline_wrap(pts, ks,2) );
    auto p0 =  bs.eval(0);
    auto p1 =  bs.eval(.2);
    auto p2 =  bs.eval(0.9);
    auto p3 =  bs.eval(1.2);
    auto p4 =  bs.eval(0.9999);

    INFO( "at 0:"<< p0 << "\n");
    INFO("at 0.2:"<< p1 << "\n");
    INFO("at 0.9:"<< p2 << "\n");
    INFO("at 1.2:"<< p3 << "\n");
    INFO("at 0.9999:"<< p4 << "\n");

    int p =  bs.degree();
    geom::ops::raise_degree(bs).swap(bs);

    REQUIRE(p == bs.degree() - 1);
    INFO("at 0:"<< bs.eval(0) << "\n");
    INFO("at 0.2:"<< bs.eval(.2) << "\n");
    INFO("at 0.9:"<< bs.eval(0.9) << "\n");
    INFO("at 1.2:"<< bs.eval(1.2) << "\n");
    INFO("at 0.9999:"<< bs.eval(0.9999) << "\n");

    REQUIRE(p0 == bs.eval(0));
    REQUIRE(p1 == bs.eval(0.2));
    REQUIRE(p2 == bs.eval(0.9));
    REQUIRE(p3 == bs.eval(1.2));
    REQUIRE(p4 == bs.eval(0.9999));
}

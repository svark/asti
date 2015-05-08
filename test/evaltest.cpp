#define CATCH_CONFIG_MAIN
#include "diego/catch/catch.hpp"
#include "bspline_x_cons.hpp"
#include "arclength.hpp"

using geom::point2d_t;
using geom::bspline;


TEST_CASE( "bspline eval test", "[bspline]" ) {

    typedef   decltype(geom::mk_stdvec(geom::point2d_t())) cpts2d_t;

    using geom::make_pt;
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
    SECTION("eval at param") {
        REQUIRE( bs.eval(0) == point2d_t(0.0));
        REQUIRE( bs.eval(.2) ==  point2d_t(0.0));
        REQUIRE( bs.eval(0.9) == point2d_t(0.0));
    }
    SECTION("eval blossom") {
        double kb[]={.2,.2,.2};
        REQUIRE( bs.blossom_eval((const double*)kb) == point2d_t(0.0));
    }

    SECTION("eval derivatives") {
        REQUIRE( geom::make_pt((bs.eval(.19) +
                                0.01 * bs.eval_derivative(1,.19) + 0.01 *0.005 *
                                bs.eval_derivative(2,0.19) ) - bs.eval(0.20)) == point2d_t(0.0));
    }

}

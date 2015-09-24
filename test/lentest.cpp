#include "point_fwd.hpp"
#include "bspline_x_cons.hpp"
#include "arclength.hpp"
#include "catch.hpp"
#include "point.hpp"
#include "testutils.hpp"
#include "conic.hpp"
#include "circle.hpp"
using geom::point2d_t;
using geom::bspline;

TEST_CASE( "bspline_arc_length", "[bspline][arclength]" ) {

    double pts[3][2] =  {{0.01, 0}, {0.0,1}, { - 1, 0}};
    double ks[] =  {0,0, 0,  1, 1, 1};
    point2d_t cpts[3];
    for(int i = 0;i < 3; ++i) {
        for(int j = 0;j < 2; ++j) {
            cpts[i][j] =  pts[i][j];
        }
    }
    bspline <point2d_t > spl2d (
        geom::make_bspline_arr(cpts, cpts + 3,
                               ks, ks + sizeof(ks) / sizeof(double), 2));
    SECTION("arc length test") {
        REQUIRE(arclength(spl2d, 0, 1) == Approx(1.5457504021));
    }
}

TEST_CASE( "circle_arc_length", "[circle][arclength]" ) {

    double pts[3][3] =  {{1, 0,1}, {0.0,1,1}, { - 1, 0,1}};

    point3d_t ps[3];
    for(int i = 0;i < 3; ++i) {
        for(int j = 0;j < 3; ++j) {
            ps[i][j] =  pts[i][j];
        }
    }
    auto const &rspl3d = 
        ( geom::make_rbspline_from_conic(geom::make_circular_arc(ps)));

    SECTION("arc length test") {
        REQUIRE(arclength(rspl3d, 0, 1) == Approx(3.14157265358));
    }
}

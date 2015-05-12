#define CATCH_CONFIG_MAIN
#include "point_fwd.hpp"
#include "bspline_x_cons.hpp"
#include "arclength.hpp"
#include "diego/catch/catch.hpp"
#include "point.hpp"
using geom::point2d_t;
using geom::bspline;
// We derive a fixture named IntegerFunctionTest from the QuickTest
// fixture.  All tests using this fixture will be automatically
// required to be quick.
TEST_CASE( "bspline arc length", "[bspline]" ) {

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

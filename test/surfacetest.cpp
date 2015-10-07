#include <vector>
#include "catch.hpp"
#include "extruded_surf.hpp"
#include "revolved_surf.hpp"
#include "circle.hpp"
#include "bspline_queries.hpp"
#include "ruled_surf.hpp"
#define _USE_MATH_DEFINES
#include <math.h>
#include "testutils.hpp"
#include "line.hpp"


TEST_CASE("surfaces","[bspline][surface]")
{
    using geom::make_pt;
    using namespace geom::qry;
    const double sqrt2 = M_SQRT1_2;
    auto circ1        = geom::make_circle(make_pt(1,0,0),
                                          make_pt(sqrt2,sqrt2,0),make_pt(0,1,0));

    auto const & bsc1 = geom::make_rbspline_from_circle(circ1);

    auto circ2       = geom::make_circle(make_pt(1,0,1),
                                         make_pt(sqrt2,sqrt2,1),make_pt(0,1,1));

    auto const & bsc2 = geom::make_rbspline_from_circle(circ2);

    SECTION("ruled surf") {
        auto const& c = make_ruled_bspline_surface(bsc1,bsc2);
        REQUIRE(num_cpts(c) == 2*num_cpts(bsc1));
        INFO("p:" <<c.eval(0,0));
    }
    SECTION("extruded surf") {
        auto const& c = make_bspline_surface(bsc1,make_vec(0,0,1));
        INFO("p:" <<c.eval(0,0));
    }
    SECTION("revolve surf") {
        auto const & l = make_line( make_pt(-3,0,0), make_vec(0,1,0));
        auto const& c = make_revolved_bspline_surf(bsc1,l, M_PI/2);
        INFO("p:" << c.eval(0,0));
    }

    return;
}

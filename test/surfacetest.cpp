#define _USE_MATH_DEFINES
#include <math.h>
#include "point.hpp"
#include <vector>
#include "catch.hpp"
#include "testutils.hpp"
#include "extruded_surf.hpp"
#include "revolved_surf.hpp"
#include "circle.hpp"
#include "ruled_surf.hpp"

TEST_CASE("surfaces","[bspline][surface]")
{
    using geom::make_pt;
    const double sqrt2 = M_SQRT1_2;
    auto circ1        = geom::make_circle(make_pt(1,0,0),
                                         make_pt(sqrt2,sqrt2,0),make_pt(0,1,0));
    auto const & bsc1 = geom::to_rational(circ1);

    auto circ2       = geom::make_circle(make_pt(1,0,1),
                                         make_pt(sqrt2,sqrt2,1),make_pt(0,1,1));
    auto const & bsc2 = geom::to_rational(circ2);

    SECTION("ruled surf") {
        auto const& c = make_ruled_bspline_surface(bsc1,bsc2);

        INFO("p:" <<c.eval(0,0));
    }
     SECTION("extruded surf") {
         auto const& c = make_bspline_surface(bsc1,make_vec(0,0,1));
         INFO("p:" <<c.eval(0,0));
     }

    return;
}


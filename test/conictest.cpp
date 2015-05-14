#include "conic.hpp"
#include "circle.hpp"
#include <vector>
#include "catch.hpp"
#include "testutils.hpp"

using geom::point2d_t;
using geom::point3d_t;
using geom::vector2d_t;

using geom::bspline;
typedef   decltype(geom::mk_stdvec(point2d_t())) cpts2d_t;
typedef   decltype(geom::mk_stdvec(point3d_t())) cpts3d_t;
using geom::make_pt;
using geom::make_vec;
TEST_CASE("conic", "[conic][2d]") {
    point2d_t pts[3];
    pts[0] = make_pt(0.0,0.0);
    pts[1] = make_pt(0.4,0.3);
    pts[2] = make_pt(1.0,0.0);
    vector2d_t v[2];
    v[0] = geom::make_vec(0,1);
    v[1] = geom::make_vec(0.4,0.5);

    geom::conic_arc<geom::point2d_t> c = geom::make_conic_arc( pts, v);

    INFO( c.eval(0) << "," << c.eval(0.5) << c.eval(1.0)<< "\n" );
    REQUIRE( c.eval(0) ==  make_pt(0, 0));
    REQUIRE( c.eval(1) ==  make_pt(1, 0));
    REQUIRE( c.eval(0.5) ==  make_pt(0.691566265,0.478915663));
    REQUIRE(c.type()== geom::ellipse);

    auto rbs ( geom::make_rbspline_from_conic < point2d_t > (c));
    INFO( rbs.eval(0) << "," << rbs.eval(0.5) << rbs.eval(1.0)<<
          rbs.eval_derivative(1, 0.5) <<
          "\n");
    REQUIRE(rbs.eval(0.5) == geom::make_pt(0.691566265,2.76626506));
    REQUIRE(rbs.eval_derivative(1, 0.5) == geom::make_vec(2.76626506,-5.18122785));
}

TEST_CASE("circle", "[circle][2d]"){
    geom::circle<geom::point2d_t> c( geom::make_pt(0.0,0.0),
                                     geom::make_pt(0,1.0));
    REQUIRE(c.eval(0) == make_pt(0, 1));
    REQUIRE(c.eval(M_PI/2) == make_pt( - 1, 0));
    REQUIRE(c.eval(M_PI) == make_pt(0,- 1));
    REQUIRE(c.eval(3*M_PI/2) ==  make_pt(1, 0));
    REQUIRE(c.eval(2 * M_PI) ==  make_pt(0, 1));
}

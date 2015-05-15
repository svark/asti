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
TEST_CASE("ellipse", "[conic][2d][ellipse]") {
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
    REQUIRE(rbs.eval(0.5) == geom::make_pt(0.691566265,0.478915663));
    REQUIRE(rbs.eval_derivative(1, 0.5) == geom::make_vec(1.85880013,0));
}

TEST_CASE("hyperbola", "[conic][2d][hyperbola]") {
    point2d_t pts[3];
    pts[0] = make_pt(0.0,0.0);
    pts[1] = make_pt(0.4,-0.3);
    pts[2] = make_pt(1.0,0.0);
    vector2d_t v[2];
    v[0] = geom::make_vec(-0.4,0.4);
    v[1] = geom::make_vec(0.4,-0.4);

    geom::conic_arc<geom::point2d_t> c = geom::make_conic_arc( pts, v);

    INFO( c.eval(0) << "," << c.eval(0.5) << c.eval(1.0)<< "\n" );
    REQUIRE( c.eval(0) ==  make_pt(0, 0));
    REQUIRE( c.eval(1) ==  make_pt(1, 0));
	auto s = make_pt(1.0237873,-0.825972281);
    REQUIRE( c.eval(0.5) ==  s);
    REQUIRE(c.type()== geom::hyperbola);

    auto rbs ( geom::make_rbspline_from_conic < point2d_t > (c));
    INFO( rbs.eval(0) << "," << rbs.eval(0.5) << rbs.eval(1.0)<<
          rbs.eval_derivative(1, 0.5) <<
          "\n");
    REQUIRE(rbs.eval(0.5) == s);
    REQUIRE(rbs.eval_derivative(1, 0.5) == geom::make_vec(0.791260076,0));
}

TEST_CASE("circle", "[circle][2d]"){
    geom::circle<geom::point2d_t> c( geom::make_pt(0.0,0.0),
                                     geom::make_pt(0,1.0));
    REQUIRE(c.eval(0) == make_pt(0, 1));
    REQUIRE(c.eval(M_PI/2) == make_pt( - 1, 0));
    REQUIRE(c.eval(M_PI) == make_pt(0,- 1));
    REQUIRE(c.eval(3*M_PI/2) ==  make_pt(1, 0));
    REQUIRE(c.eval(2 * M_PI) ==  make_pt(0, 1));

	auto ac = geom::make_circle( geom::make_pt(0,1.0), geom::make_pt(1,0), geom::make_pt(-1,0));
	REQUIRE( c.getCenter() == ac.getCenter());
	REQUIRE( c.getRadius() == Approx(ac.getRadius()));
	REQUIRE( len (c.getPlaneNormal()  + ac.getPlaneNormal() ) == Approx(0));

	auto const & rbsc = geom::to_rational( c );
	REQUIRE(rbsc.eval(0 ) == make_pt(0, 1));
    REQUIRE(rbsc.eval(2*M_PI/3.0 ) == make_pt(-0.866025404,-0.5));
    REQUIRE(rbsc.eval(2*M_PI/6.0) == make_pt(-0.866025404,0.5));
    REQUIRE(rbsc.eval(2*M_PI/2) ==  make_pt(0, -1));
    REQUIRE(rbsc.eval(4.0*M_PI/3) ==  make_pt(0.866025404,-0.5));
}

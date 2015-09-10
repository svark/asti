#include "point.hpp"
#include <vector>
#include "catch.hpp"
#include "bspline.hpp"
#include "testutils.hpp"
#include "periodic_bspline_cons.hpp"
#include "reduce_degree.hpp"

TEST_CASE("reduce_degree", "[periodic_bspline][bspline][degree][approximation]")
{
    cpts2d_t pts(5);
    pts[0] = make_pt(0.0,0.0);
    pts[1] = make_pt(0.4,0.3);
    pts[2] = make_pt(0.2,0.8);
    pts[3] = make_pt(0.1, 0.4);
    pts[4] = make_pt(-0.2, 0.4);

    //pts[4] = make_pt(0.0,0.0);
    std::vector<double> ks(6);
    ks[0] = 0.0; ks[1] = 0.3; ks[2] = 0.6;
    ks[3] = 0.8; ks[4] = 0.9; ks[5] = 1.0;
    geom::periodic_bspline<geom::point2d_t> bs( geom::make_periodic_bspline_wrap(pts, ks,3) );
    auto const  & nc = geom::ops::reduce_degree( bs, 2);
	REQUIRE( geom::qry::is_periodic(geom::qry::get_spline(nc)) );
    auto p0 =  bs.eval(0);
    auto p1 =  bs.eval(.2);
    auto p2 =  bs.eval(0.9);
    auto p3 =  bs.eval(1.2);
    auto p4 =  bs.eval(0.9999999);

   
    auto q0 =  nc.eval(0);
    auto q1 =  nc.eval(.2);
    auto q2 =  nc.eval(0.9);
    auto q3 =  nc.eval(1.2);
    auto q4 =  nc.eval(0.9999999);
	double dist = len(p0-q0);
	REQUIRE(dist== Approx(0));
	dist += len(p1-q1);
	dist += len(p2-q2);
	dist += len(p3-q3);
	double diste = len(p4-q4);
	REQUIRE(diste == Approx(0.0));
	dist += diste;
	REQUIRE(dist < 0.18);

}

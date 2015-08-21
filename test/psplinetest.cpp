// #define CATCH_CONFIG_MAIN
#include "point.hpp"
#include <vector>
#include "periodic_bspline_cons.hpp"
#include "bspline_queries.hpp"
#include "bspline_ops.hpp"
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
TEST_CASE("periodic_bspline_basic", "[periodic][2d][eval][derivative][blossom][bspline]" ) {

    cpts2d_t pts(4);

    pts[0] = make_pt(0.0,0.0);
    pts[1] = make_pt(0.4,0.3);
    pts[2] = make_pt(0.2,0.8);
    pts[3] = make_pt(-0.2, 0.4);
    //pts[4] = make_pt(0.0,0.0);
    std::vector<double> ks(5);
    ks[0] = 0.0; ks[1] = 0.3; ks[2] = 0.6;
    ks[3] = 0.8; ks[4] = 1.0;
    geom::periodic_bspline<point2d_t> bs =
        geom::make_periodic_bspline_wrap(pts, ks, 2);

    INFO( "at 0:"<< bs.eval(0) << "\n");
    REQUIRE(bs.eval(0.0) == make_pt(0.04,0.64));
    INFO("at 0.2:"<< bs.eval(.2) << "\n");
    INFO("at 0.9:"<< bs.eval(0.9) << "\n");
    REQUIRE(bs.eval(0.2) == make_pt(-0.12888889, 0.33777778));
    REQUIRE(bs.eval(0.9) == make_pt(0.185,0.6975));
    INFO("at 1.2:"<< bs.eval(1.2) << "\n");
    INFO("at 0.9999:"<< bs.eval(.9999) << "\n");
    REQUIRE(bs.periodic_param(1.2) == Approx(0.2));
    REQUIRE(bs.periodic_param(-0.2) == Approx(0.8));
    REQUIRE(bs.periodic_param(-1.2) == Approx(0.8));
    REQUIRE(geom::qry::is_periodic(bs.spline()));
    REQUIRE(len(bs.eval_derivative(1, 0.9 - 1e-7) - bs.eval_derivative(1, 0.9)) <
            1e-7 * len(bs.eval_derivative(2, 0.9))
        );


}

TEST_CASE("periodic_bspline_rational", "[reparametrise][split][rational][periodic][2d][3d]"){
    cpts3d_t pts(4);
    pts[0] = make_pt(0.0,0.0,1);
    pts[1] = make_pt(0.4,0.3,1);
    pts[2] = make_pt(0.2,0.8,1);
    pts[3] = make_pt(-0.2, 0.4,1);
    std::vector<double> ks(5);
    ks[0] = 0.0; ks[1] = 0.3; ks[2] = 0.6;
    ks[3] = 0.8; ks[4] = 1.0;
    auto  rbs = geom::make_rbspline(std::move(geom::make_periodic_bspline_wrap(pts, ks,2)));
    auto const & pc = geom::ops::rotate_base_knot(rbs.spline(), 3);
    auto rpc = geom::ops::reparametrize(pc,0,1);
    auto rbs1 = geom::ops::split_periodic_curve(rbs, 0.1);
    point2d_t pt01 ( rbs1.eval(0.1) );
    REQUIRE(pt01 == make_pt(-0.082222222,0.48444444));
    auto const &vs = rbs1.eval_derivatives(2,0.1);
    vector2d_t v1 =  vs.front();
    vector2d_t v2 =  vs.back();
    REQUIRE(v1 == make_vec(pt01));
    REQUIRE(v2 == make_vec(7.55555556,0.888888889));

}

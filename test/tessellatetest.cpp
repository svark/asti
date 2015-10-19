#include "point.hpp"
#include <vector>
#include "catch.hpp"
#include "testutils.hpp"
#include "tessellate.hpp"
#include <functional>
#include "tol.hpp"
#include "bspline_queries.hpp"
#include "rational_bspline_cons.hpp"
#include "circle.hpp"
#include "line.hpp"

template <class SplineType>
void
tess_is_close(const SplineType& bs, const double tess_tol)
{
    using geom::qry::num_cpts;
	using geom::qry::start_param;
	using geom::qry::end_param;

    geom::polyline<geom::point2d_t> const &pl = geom::tessellate(bs, tess_tol);
	auto s = start_param(bs);
    auto w = end_param(bs) - s;
    size_t d = num_cpts(pl)-1;
    auto f = w/d;
    point2d_t pp;
    for(size_t i = 0 ; i < pl.control_points().size() ; ++i)
    {
        auto p  =  pl.control_points()[i];

        REQUIRE(tol::small(len(p - bs.eval(s + i*f))));

        if( i != 0 ) {
            REQUIRE(geom::distance_from_point_to_line_seg(
                        geom::make_line_seg(p,pp) , bs.eval(s + (i-0.5)*f)) < tess_tol);
        }
        pp = p;
    }
}

TEST_CASE("tessellate bspline",  "[bspline][linear][approximation]") {


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
    geom::bspline<geom::point2d_t> bs((pts),
                                      (ks),2);

    const double tess_tol = .2;

    SECTION("bspline")
    {
        tess_is_close(bs, tess_tol);
    }
    SECTION("rational_bspline")
    {
        geom::polyline<geom::point2d_t> const &pl = geom::tessellate(bs, tess_tol);
        std::vector<double> weights(5,1.0);
        auto const &rbs = geom::make_rbspline(std::move(pts), std::move(weights), std::move(ks),  2);

        geom::polyline<geom::point2d_t> const &pl2 = geom::tessellate(rbs, tess_tol);
        REQUIRE(pl2 == pl);

    }
    SECTION("circle")
    {
        geom::point2d_t p1 = make_pt(-1,1.5),p2 = make_pt(0,1),p3 = make_pt(1,0);
        auto const &rbs = geom::make_rbspline_from_circle(geom::make_circle( p1, p2, p3 ));
        tess_is_close(rbs, 0.1);
    }
}

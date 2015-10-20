#include "point.hpp"
#include <vector>
#include "catch.hpp"
#include "bspline.hpp"
#include "testutils.hpp"
#include "split_into_bezier_patches.hpp"
#include "smat.hpp"
#include "bspline_queries.hpp"
#include "conic.hpp"
#include "circle.hpp"

TEST_CASE("bezier_test","[bezier][bspline][split]") {

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

    auto spls = geom::ops::split_into_bezier_patches(bs);
    double us[] = {0.1,0.6,0.9};
    int i = 0;
    for(auto s : spls)
    {
        auto const & p1 =  s.eval(us[i]);
        auto const & p2 =  bs.eval(us[i]);
        INFO( "bez:" << p1 << "spl:" << p2 << "\n" );
        REQUIRE( p1 == p2);
        ++i;
    }
}

TEST_CASE("bezier_test1", "[bezier][rational bspline][split][circle]" ){
	using geom::make_pt;
    geom::circle<geom::point2d_t> c( make_pt(0.0,0.0),
                                     make_pt(0,2.0));

	
    auto const & rbsc = geom::make_rbspline_from_circle( c );
   
    auto const &bps = geom::ops::split_into_bezier_patches(rbsc);
    double v = 0;
    for(auto bp : bps)
    {
        double u = (0.25*geom::qry::start_param(bp) + 0.75*geom::qry::end_param(bp));
        INFO("u:" << u);
        REQUIRE( bp.eval(u) == rbsc.eval(u));
        v = u;
    }
}
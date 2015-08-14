//-*- mode:c++ -*-
#include "point.hpp"
#include <vector>
#include "catch.hpp"
#include "testutils.hpp"
#include "trim_extend_join.hpp"
#include "subdivide_curve.hpp"
#include "reverse_curve.hpp"
#include "spline_approx.hpp"
#include "reparametrize.hpp"

TEST_CASE("trimextendjointest",  "[bspline][trim][extend][join]")
{
    cpts2d_t pts(5);
    pts[0] = make_pt(0.0,0.0);
    pts[1] = make_pt(1.0,1.0);
    pts[2] = make_pt(1.5,0.3);
    pts[3] = make_pt(1.8,0.1);
    pts[4] = make_pt(2.0,0.0);
    std::vector<double> ks(3+3+2);
    ks[0] = ks[1] = ks[2] = 0.0;
    ks[3] = 0.5; ks[4] = 0.8;
    ks[5] = ks[6] = ks[7] = 1.2;
    geom::bspline<geom::point2d_t> bs(std::move(pts), std::move(ks),2);

    SECTION("trim")
    {
        auto c1 ( geom::ops::trim_curve(bs,0, 0.6) ) ;
        auto c2 ( geom::ops::trim_curve(bs,0.6,1.2) );
        for(int i =0;i < 6;i+=2)
        {
            auto const &p =  c1.eval(i/10.0);
            auto fp =  geom::ops::foot_param(bs,p);
            REQUIRE(len(bs.eval(fp) - p) == Approx(0.0));
        }
        for(int i =6;i < 12;i+=2)
        {
            auto const &p =  c2.eval(i/10.0);
            auto fp =  geom::ops::foot_param(bs,p);
            REQUIRE(len(bs.eval(fp) - p) == Approx(0.0));
        }
    }
    SECTION("extend")
    {
        auto c1 ( geom::ops::trim_curve(bs,0.1, 1.0) ) ;
        auto c3 = geom::ops::extend_curve_start(geom::ops::extend_curve_end(c1,0.1),0.1);
         for(int i =0;i < 6;i+=2)
        {
            auto const &p =  c3.eval(i/10.0);
            auto fp =  geom::ops::foot_param(bs,p);
            REQUIRE(len(bs.eval(fp) - p) == Approx(0.0));
        }
        for(int i =7;i < 12;i+=2)
        {
            auto const &p =  c3.eval(i/10.0);
            auto fp =  geom::ops::foot_param(bs,p);
            REQUIRE(len(bs.eval(fp) - p) == Approx(0.0));
        }
    }
    SECTION("join")
    {
        auto c1 ( geom::ops::trim_curve(bs,0, 0.6) ) ;
        auto c2 ( geom::ops::trim_curve(bs,0.601,1.2) );
        auto c3 ( geom::ops::reparametrize( geom::ops::join_starts(
                                                geom::ops::reverse_curve(c1), c2, 1) ,0,1.2)) ;
        for(int i =0;i < 6;i+=2)
        {
            auto const &p =  c3.eval(i/10.0);
            auto fp =  geom::ops::foot_param(bs,p);
            REQUIRE(len(bs.eval(fp) - p) == Approx(0.0));
        }
        for(int i =7;i < 12;i+=2)
        {
            auto const &p =  c3.eval(i/10.0);
            auto fp =  geom::ops::foot_param(bs,p);
            REQUIRE(len(bs.eval(fp) - p) == Approx(0.0));
        }

        auto c4( geom::ops::reparametrize(
                     geom::ops::extend_curve_end_to_pt(
                         geom::ops::reverse_curve(c1), bs.eval(0.601)), 0,.601));
        for(int i =0;i < 6;i+=2)
        {
            auto const &p =  c4.eval(i/10.0);
            auto fp =  geom::ops::foot_param(bs,p);
            REQUIRE(len(bs.eval(fp) - p) == Approx(0.0));
        }
        REQUIRE(len(bs.eval(0.601) - c4.eval(0.601)) == Approx(0.0));

    }
}

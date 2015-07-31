#include "point.hpp"
#include <vector>
#include "catch.hpp"
#include "testutils.hpp"
#include "spline_approx.hpp"
#include "spline_interp.hpp"
#include <functional>

TEST_CASE("approxtestbasic",  "[bspline][approximation]") {

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

  const double ts[] = {-1,-1,-1,-1,0.5,1,1,1,1};
  auto cube = [](double x) { return x*x*x -3*x+1;};
  auto spl0 = geom::ops::cubic_approx1d(
      std::function<double(double)>(cube),
      std::vector<double>(ts, ts + sizeof(ts) / sizeof(ts[0])));
  INFO( "cube approx at -1,0.1,0,1:"
            << spl0.eval(-1) << ","
            << spl0.eval(0.1) << ',' << spl0.eval(0) << ','
        << spl0.eval(1) << "\n") ;

  REQUIRE(spl0.eval( - 1) == Approx(3));
  REQUIRE(spl0.eval( 0.1) == Approx(0.701));
  REQUIRE(spl0.eval( 0)   == Approx(1));
  REQUIRE(spl0.eval( 1)   == Approx( - 1));

  auto p = make_pt(0.68,0.586);
  std::vector<double> newknots(bs.knots().size() + 2);
  double arr[] = {0.225,0.75};
  std::merge(bs.knots().cbegin(), bs.knots().cend(), arr, arr + 2, newknots.begin());
  newknots.insert(newknots.begin(), 0.0);
  newknots.insert(newknots.end(), 1.0);
  auto v_vdash = [&p,&bs]
    (double u) -> double
      {
    auto v =   (bs.eval(u) - p);
    auto vdash =  bs.eval_derivative(1,u);
    return dot(v, vdash);
  };

  auto spl2 = geom::ops::cubic_approx1d(
      std::function<double(double)>(v_vdash), newknots);
  auto pt_ = spl2.eval(-1);
  INFO(  "\n dist.dist'  at 0.1,0.2,1:" << v_vdash(0.25)
         << ',' << v_vdash(0.2) << ',' << v_vdash(.5) << "\n");
  INFO(  "\n dist.dist' approx at 0.1,0.2,1:"
             << spl2.eval(0.25) << ',' << spl2.eval(0.2) << ','
         << spl2.eval(.5) << "\n");

  REQUIRE( fabs(v_vdash(0.2) - spl2.eval(0.2)) < 0.05);
  REQUIRE( fabs(v_vdash(0.2) - spl2.eval(0.2)) < 0.05);
  REQUIRE( fabs(v_vdash(0.5) - spl2.eval(0.5)) < 0.1);

  auto fp =  geom::ops::foot_param(bs,p);
  INFO( "soln NEAR param .2::" << fp);
  REQUIRE(fp ==  Approx(0.199843));
}

TEST_CASE("cubic hermite interp", "[bspline][interpolation][approximation]")
{
  std::vector<double> ps(6);
  ps[0] = 0;ps[1] = 2;ps[2] = 2.5;ps[3] = 3;ps[4] = 3.25;ps[5] = 3.1;
  geom::interpolation_options_t opts;
  opts.parametrization = geom::centripetal_length;
  opts.end_conditions =  geom::not_a_knot;
  SECTION("not a knot") {
  auto bs= geom::piecewise_cubic_hermite_interp(ps.begin(), ps.end(),opts,
                                          std::vector<double>());
  INFO( "v:" <<  bs.eval(0.2) << ","  << bs.eval(0.3806)  << ","<< bs.eval(0.4) << "," << bs.eval(0.5) << "," <<bs.eval(0.6) << "," << bs.eval(0.9) << "," << bs.eval(1.1) << "\n" );
 //1.31061,1.99999,2.05698,2.32315,2.57359,3.24986,2.90523
  REQUIRE(bs.eval(0) == Approx(0.0));
  REQUIRE(bs.eval(0.2) == Approx(1.31061));
  REQUIRE(bs.eval(0.3806) == Approx(2));
  REQUIRE(bs.eval(0.4) == Approx(2.05698));
  REQUIRE(bs.eval(0.5) == Approx(2.32315));
  REQUIRE(bs.eval(0.6) == Approx(2.57359));
  REQUIRE(bs.eval(0.9) == Approx(3.24986));
  REQUIRE(bs.eval(1.1) == Approx(2.90523));

  }
  opts.end_conditions =  geom::vanishing_double_derivatives;
  SECTION("vanishing double derivatives") {
  auto bs= geom::piecewise_cubic_hermite_interp(ps.begin(), ps.end(),opts,
                                           std::vector<double>());
  INFO( "v:" <<  bs.eval(0.2) << "," << bs.eval(0.4) << "," << bs.eval(0.5) << "," <<bs.eval(0.6) << "," << bs.eval(0.9) << "," << bs.eval(1.1) << "\n" );
 // 0.802717,2.08463,2.36436,2.56361,3.24958,2.95042
  REQUIRE(bs.eval(0) == Approx(0.0));
  REQUIRE(bs.eval(0.2) == Approx(0.802717));
  REQUIRE(bs.eval(0.3806) == Approx(2));
  REQUIRE(bs.eval(0.4) == Approx(2.08463));
  REQUIRE(bs.eval(0.5) == Approx(2.36436));
  REQUIRE(bs.eval(0.6) == Approx(2.56361));
//  REQUIRE(bs.eval(0.8957) == Approx( 3.2499998098 ));
  REQUIRE(bs.eval(0.9) == Approx(3.24958));
  REQUIRE(bs.eval(1.1) == Approx(2.95042));

  }
   opts.end_conditions =  geom::parabolic_blending;
   SECTION("parabolic_blending") {
  auto bs= geom::piecewise_cubic_hermite_interp(ps.begin(), ps.end(),opts,
                                          std::vector<double>());
  INFO( "v:" << bs.eval(0.2) << "," << bs.eval(0.4) << "," << bs.eval(0.5) << "," <<bs.eval(0.6) << "," << bs.eval(0.9) << "," << bs.eval(1.1) << "\n" );

  REQUIRE(bs.eval(0) == Approx(0.0));
  REQUIRE(bs.eval(0.2) == Approx( 0.997261));
  REQUIRE(bs.eval(0.3806) == Approx(2));
  REQUIRE(bs.eval(0.4) == Approx(2.07433));
  REQUIRE(bs.eval(0.5) == Approx( 2.3501));
  REQUIRE(bs.eval(0.6) == Approx(2.5658));
  REQUIRE(bs.eval(0.9) == Approx(3.24719));
  REQUIRE(bs.eval(1.1) == Approx(3.54124 ));
  //REQUIRE(bs.eval(0.8957) == Approx(3.2499998098));
  }
   //ps.back() = 3.0;
   ps.push_back(2.5);
   ps.push_back(2.0);
   ps.push_back(0.0);
   SECTION("periodic") {
       opts.end_conditions=geom::periodic;
        auto bs= geom::piecewise_cubic_hermite_interp_periodic(
            ps.begin(), ps.end(),opts, std::vector<double>());
  INFO( "v:" << bs.eval(0) << "," << bs.eval(0.2) << "," << bs.eval(0.3806) << "," << bs.eval(0.4) << "," << bs.eval(0.5) << "," <<bs.eval(0.6) << "," << bs.eval(0.9) << "," << bs.eval(1.1) << "\n" );
	//1.43328,0.000132328,2.32469,2.40941,2.96996,3.23082,2.20675,0.218459
  REQUIRE(bs.eval(0) == Approx(0.0));
  REQUIRE(bs.eval(0.2) == Approx(1.87719));
  REQUIRE(bs.eval(0.3806) == Approx(2.76392));
  REQUIRE(bs.eval(0.4) == Approx(2.86035));
  REQUIRE(bs.eval(0.5) == Approx(3.2487 ));
  REQUIRE(bs.eval(0.6) == Approx(2.90866));
  //REQUIRE(bs.eval(0.8957) == Approx(0));
  REQUIRE(bs.eval(0.9) == Approx(0.705953 ));
  REQUIRE(bs.eval(1.1) == Approx(0.699618));
  REQUIRE( geom::ops::is_periodic(bs.spline()));
  //INFO ("is periodic" << std::boolalpha << ":" << geom::ops::is_periodic(bs.spline()) << "\n\n");
  for( auto k : bs.knots() )
  {
      INFO ("at "<< k << ": " << bs.eval(k) << "\n");
  }
  INFO ("------------\n");
   REQUIRE(bs.eval(0) == Approx(0.0));
   }
}

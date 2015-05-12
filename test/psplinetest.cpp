#define CATCH_CONFIG_MAIN
#include "point.hpp"
#include <vector>
#include "periodic_bspline_cons.hpp"
#include <ostream>
#include "bspline_queries.hpp"
namespace std {
ostream &
operator<<(ostream& os,
           const ::geom::pt_t<2>& pt)
{
  os << "(" <<pt.p[0] <<"," << pt.p[1] << ")";
  return os;
}
ostream &
operator<<(ostream& os,
           const ::geom::vec_t<2>& v)
{
  os << "(" <<v.v[0] <<"," << v.v[1] << ")";
  return os;
}
}
#include "diego/catch/catch.hpp"
using geom::point2d_t;
using geom::bspline;

TEST_CASE( "periodic_bspline eval test", "[periodic_bspline][2d][eval][derivative][blossom]" ) {
    typedef   decltype(geom::mk_stdvec(geom::point2d_t())) cpts2d_t;

  cpts2d_t pts(4);
  using geom::make_pt;
  pts[0] = make_pt(0.0,0.0);
  pts[1] = make_pt(0.4,0.3);
  pts[2] = make_pt(0.2,0.8);
  pts[3] = make_pt(-0.2, 0.4);
  //pts[4] = make_pt(0.0,0.0);
  std::vector<double> ks(5);
  ks[0] = 0.0; ks[1] = 0.3; ks[2] = 0.6;
  ks[3] = 0.8; ks[4] = 1.0;
  geom::periodic_bspline<geom::point2d_t> bs = geom::make_periodic_bspline_wrap(pts, ks, 2);
  
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
  bs.optimize();
  INFO( "optimized\nat 0:"<< bs.eval(0) << "\n");
  REQUIRE(bs.eval(0.0) == make_pt(0.04,0.64));
  REQUIRE(bs.eval(0.2) == make_pt(-0.128888889,0.337777778));
  REQUIRE(bs.eval(0.9) == make_pt(0.185,0.6975));
  REQUIRE(geom::ops::is_periodic(bs.spline()));
  REQUIRE(len(bs.eval_derivative(1, 1 - 1e-7) - bs.eval_derivative(1, 1.0)) < 1e-5);
}

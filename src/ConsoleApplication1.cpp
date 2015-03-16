// ConsoleApplication1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <iostream>
#include <memory>
#include <utility>
#include <vector>
#include <list>
#include <functional>
#include "point.hpp"
#include "spline_interp.hpp"
#include "circle.hpp"

std::ostream& operator<<( std::ostream& os, const geom::point2d_t& pt)
{
  os << "(" <<pt.p[0] <<"," << pt.p[1] << ")";
  return os;
}

std::ostream& operator<<( std::ostream& os, const geom::vector2d_t& pt)
{
  os << "(" <<pt.v[0] <<"," << pt.v[1] << ")";
  return os;
}
#include "geom_exception.hpp"
#include "bspline.hpp"
#include "periodic_bspline.hpp"
#include "bspline_ops.hpp"
#include "trim_extend_join.hpp"
#include "subdivide_curve.hpp"
#include "reverse_curve.hpp"
#include "reparametrize.hpp"
#include "bspline_queries.hpp"
#include "rational_bspline.hpp"
#include "conic.hpp"
#include "rotate_base_knot.hpp"
double cube(double u) { double x = u ; return x*x*x -3*x+1;}

int _tmain(int argc, _TCHAR* argv[])
{
  typedef   geom::point_iter_traits<geom::point2d_t*>::PointContT cpts2d_t;
  typedef   geom::point_iter_traits<geom::point3d_t*>::PointContT cpts3d_t;
  using geom::make_pt;
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
  std::cout << "at 0:"<< bs.eval(0) << "\n";
  std::cout <<"at 0.2"<< bs.eval(.2) << "\n";
  std::cout <<"at 0.9"<< bs.eval(0.9) << "\n";
  double kb[]={.2,.2,.2};
  std::cout <<"at 0.2"<< bs.blossom_eval((const double*)kb) << "\n";
  std::cout << "diff at .22" << geom::make_pt((bs.eval(.19) +
                                               0.01 * bs.eval_derivative(1,.19) + 0.01 *0.005 *
                                               bs.eval_derivative(2,0.19) ) - bs.eval(0.20)) <<"\n";
  auto p = make_pt(0.68,0.586);

  auto fn = [&p,&bs]
    (double u) -> std::tuple<double,double,double,double> {
    auto v         = (bs.eval(u) - p);
    auto vdash     = bs.eval_derivative(1,u);
    auto vdashdash = bs.eval_derivative(2,u);
    return std::make_tuple(dot(v, vdash),
                           dot(v, vdashdash ) + dot(vdash,vdash),
                           3 * dot(vdash, vdashdash),
                           3 * dot(vdashdash,vdashdash)
                           );
  };

  auto v_vdash = [&p,&bs]
    (double u) -> double {
    auto v =   (bs.eval(u) - p);
    auto vdash =  bs.eval_derivative(1,u);
    return dot(v, vdash);
  };

  double f0,f1,f2,f3;
  std::tie(f0,f1,f2,f3) = fn(0.19);

  std::cout << '(' << f0 << ',' << f1  << ')' << '\n';
  std::cout << "diff at .20, 0.19:"
            << f0 + 0.01 * f1 + (0.01 * 0.01* f2) / 2.0 +
      (0.01 * 0.01 *0.01 * f3) / 6.0
      - std::get<0>(fn(0.20)) << "\n";


  std::cout << "diff at .21, 0.19"
            << geom::make_pt((bs.eval(.19) + 0.02 *
                              bs.eval_derivative(1,.19)
                              + 0.02 *0.01 *
                              bs.eval_derivative(2,0.19))
                             - bs.eval(0.21))
            << "\n" ;
  std::cout <<"dot:" << dot(geom::make_vec( make_pt(1.0 , 0)) , geom::make_vec(make_pt(0.0,1.0) )) << "\n";

  std::cout << '(' << f0 <<',' <<  f1  << ')' << '\n';
  const double ts[] = {-1,-1,-1,-1,0.5,1,1,1,1};
  auto spl0 = geom::bspline_ops::cubic_approx1d(
      std::function<double(double)>(cube),
      std::vector<double>(ts, ts + sizeof(ts) / sizeof(ts[0])));
  std::cout << "cube approx at -1,0.1,0,1:"
            << spl0.eval(-1) << ","
            << spl0.eval(0.1) << ',' << spl0.eval(0) << ','
            << spl0.eval(1) << "\n";

  std::vector<double> newknots(bs.knots().size() + 2);
  double arr[] = {0.225,0.75};
  std::merge(bs.knots().cbegin(), bs.knots().cend(), arr, arr + 2, newknots.begin());
  newknots.insert(newknots.begin(), 0.0);
  newknots.insert(newknots.end(), 1.0);
  auto spl2 = geom::bspline_ops::cubic_approx1d(
      std::function<double(double)>(v_vdash), newknots);
  try {
    auto pt_ = spl2.eval(-1);
    std::cout << "should not reach here\n";
  }catch(geom::geom_exception e) {
    std::cout << e.what() << "\n";
  }
  std:: cout << "\n dist.dist'  at 0.1,0.2,1:" << v_vdash(0.25)
             << ',' << v_vdash(0.2) << ',' << v_vdash(.5) << "\n";
  std:: cout << "\n dist.dist' approx at 0.1,0.2,1:"
             << spl2.eval(0.25) << ',' << spl2.eval(0.2) << ','
             << spl2.eval(.5) << "\n";
  std::cout << "soln NEAR param .2::" << geom::bspline_ops::foot_param(bs,p);
  auto spls = geom::bspline_ops::split_into_bezier_patches(bs);
  double us[] = {0.1,0.6,0.9};
  int i=0;
  for(auto s : spls)
  {
      std::cout << "bez:" << s.eval(us[i]) << "spl:" << bs.eval(us[i]) << "\n";
      ++i;
  }
  {
      auto spl = geom::bspline_ops::extend_curve_start(bs,0.2);
      std::cout << "at 0.15" << spl.eval(0.15) << ",actual:"<< bs.eval(0.15) <<"\n";
      std::cout << "at -0.15" << spl.eval(-0.15) << ",actual:"<< bs.eval(-0.15) <<"\n";
      std::cout << "at -0.2" << spl.eval(-0.2) << ",actual:"<< bs.eval(-0.2) <<"\n";
  }
  {
  cpts2d_t pts(4);
  pts[0] = make_pt(0.0,0.0);
  pts[1] = make_pt(0.4,0.3);
  pts[2] = make_pt(0.2,0.8);
  pts[3] = make_pt(-0.2, 0.4);
  //pts[4] = make_pt(0.0,0.0);
  std::vector<double> ks(5);
  ks[0] = 0.0; ks[1] = 0.3; ks[2] = 0.6;
  ks[3] = 0.8; ks[4] = 1.0;
  geom::periodic_bspline<geom::point2d_t> bs= geom::make_periodic_bspline_wrap(pts, ks,2);
  std::cout << "at 0:"<< bs.eval(0) << "\n";
  std::cout <<"at 0.2:"<< bs.eval(.2) << "\n";
  std::cout <<"at 0.9:"<< bs.eval(0.9) << "\n";
  std::cout <<"at 1.2:"<< bs.eval(1.2) << "\n";
  std::cout <<"at 0.9999:"<< bs.eval(.9999) << "\n";
  bs = geom::bspline_ops::raise_degree(bs);
   std::cout << "at 0:"<< bs.eval(0) << "\n";
  std::cout <<"at 0.2:"<< bs.eval(.2) << "\n";
  std::cout <<"at 0.9:"<< bs.eval(0.9) << "\n";
  std::cout <<"at 1.2:"<< bs.eval(1.2) << "\n";
  std::cout <<"at 0.9999:"<< bs.eval(0.9999) << "\n";
  bs.optimize();
  std::cout << "optimized\nat 0:"<< bs.eval(0) << "\n";
  std::cout <<"at 0.2:"<< bs.eval(.2) << "\n";
  std::cout <<"at 0.9:"<< bs.eval(0.9) << "\n";
  std::cout <<"at 1.2:"<< bs.eval(1.2) << "\n";
  std::cout <<"at 0.9999:"<< bs.eval(.9999) << "\n";
  bs = geom::bspline_ops::raise_degree(bs);
   std::cout << "at 0:"<< bs.eval(0) << "\n";
  std::cout <<"at 0.2:"<< bs.eval(.2) << "\n";
  std::cout <<"at 0.9:"<< bs.eval(0.9) << "\n";
  std::cout <<"at 1.2:"<< bs.eval(1.2) << "\n";
  std::cout <<"at 0.9999:"<< bs.eval(0.9999) << "\n";
  }

  std::vector<double> ps(6);
  ps[0] = 0;ps[1] = 2;ps[2] = 2.5;ps[3] = 3;ps[4] = 3.25;ps[5] = 3.1;
  geom::interpolation_options_t opts;
  opts.parametrization = geom::centripetal_length;
  opts.end_conditions =  geom::not_a_knot;
  {
  auto bs= geom::piecewise_cubic_hermite_interp(ps.begin(), ps.end(),opts,
                                          std::vector<double>());
  std::cout << "at 0:"<< bs.eval(0) << "\n";
  std::cout <<"at 0.2:"<< bs.eval(.2) << "\n";
  std::cout <<"at 0.3806:"<< bs.eval(.3806) << "\n";
  std::cout <<"at 0.4:"<< bs.eval(0.4) << "\n";
  std::cout <<"at 0.5:"<< bs.eval(0.5) << "\n";
  std::cout <<"at 0.6:"<< bs.eval(0.6) << "\n";
  std::cout <<"at 0.9:"<< bs.eval(0.9) << "\n";
  std::cout <<"at 1.1:"<< bs.eval(1.1) << "\n";
  std::cout <<"at 0.8957:"<< bs.eval(0.8957) << "\n";
  }
  opts.end_conditions =  geom::vanishing_double_derivatives;
   {
  auto bs= geom::piecewise_cubic_hermite_interp(ps.begin(), ps.end(),opts,
                                          std::vector<double>());
  std::cout << "at 0:"<< bs.eval(0) << "\n";
  std::cout <<"at 0.2:"<< bs.eval(.2) << "\n";
  std::cout <<"at 0.3806:"<< bs.eval(.3806) << "\n";
  std::cout <<"at 0.4:"<< bs.eval(0.4) << "\n";
  std::cout <<"at 0.5:"<< bs.eval(0.5) << "\n";
  std::cout <<"at 0.6:"<< bs.eval(0.6) << "\n";
  std::cout <<"at 0.9:"<< bs.eval(0.9) << "\n";
  std::cout <<"at 1.1:"<< bs.eval(1.1) << "\n";
  std::cout <<"at 0.8957:"<< bs.eval(0.8957) << "\n";
  }
   opts.end_conditions =  geom::parabolic_blending;
   {
  auto bs= geom::piecewise_cubic_hermite_interp(ps.begin(), ps.end(),opts,
                                          std::vector<double>());
  std::cout << "at 0:"<< bs.eval(0) << "\n";
  std::cout <<"at 0.2:"<< bs.eval(.2) << "\n";
  std::cout <<"at 0.3806:"<< bs.eval(.3806) << "\n";
  std::cout <<"at 0.4:"<< bs.eval(0.4) << "\n";
  std::cout <<"at 0.5:"<< bs.eval(0.5) << "\n";
  std::cout <<"at 0.6:"<< bs.eval(0.6) << "\n";
  std::cout <<"at 0.9:"<< bs.eval(0.9) << "\n";
  std::cout <<"at 1.1:"<< bs.eval(1.1) << "\n";
  std::cout <<"at 0.8957:"<< bs.eval(0.8957) << "\n";
  }
   ps.push_back(2.5);
   ps.push_back(2.0);
   ps.push_back(0.0);
   {
       opts.end_conditions=geom::periodic;
        auto bs= geom::piecewise_cubic_hermite_interp_periodic(ps.begin(), ps.end(),opts,
                                          std::vector<double>());
  std::cout << "at 0:"<< bs.eval(0) << "\n";
  std::cout <<"at 0.2:"<< bs.eval(.2) << "\n";
  std::cout <<"at 0.3806:"<< bs.eval(.3806) << "\n";
  std::cout <<"at 0.4:"<< bs.eval(0.4) << "\n";
  std::cout <<"at 0.5:"<< bs.eval(0.5) << "\n";
  std::cout <<"at 0.6:"<< bs.eval(0.6) << "\n";
  std::cout <<"at 0.9:"<< bs.eval(0.9) << "\n";
  std::cout <<"at 1.1:"<< bs.eval(1.1) << "\n";
  std::cout <<"at 0.8957:"<< bs.eval(0.8957) << "\n";
  std::cout <<"is periodic" << std::boolalpha << ":" << geom::bspline_ops::is_periodic(bs.spline()) << "\n\n";
  for( auto k : bs.knots() ) 
  {
      std::cout <<"at "<< k << ": " << bs.eval(k) << "\n";
  }
  std::cout <<"------------\n";
   }
  //auto spls2 = geom::bspline_ops::split_into_bezier_patches_hard(bs);
  //bool ok = std::equal(spls.begin(), spls.end(),spls2.begin() );
  //std::cout << std::boolalpha << ok;
   {
       auto &p = geom::bspline_ops::split_open_curve(bs,0.3);
       auto &p_ = geom::bspline_ops::split_open_curve(geom::bspline_ops::reverse_curve(bs),-0.3);
       geom::bspline<geom::point2d_t> s1 = geom::bspline_ops::reverse_curve(p.first),s2 = p.second;
       auto const &p2 = geom::bspline_ops::reparametrize_start(
           geom::bspline_ops::join_starts(s1,s2,1), 0);
       

       std::cout << "----------\nat 0:"<< p.first.eval(0) << "\n";
  std::cout <<"at 0.2:"<< p.first.eval(.2) << "\n";
  std::cout <<"at 0.3:"<< p.first.eval(.3) << "\n";
  std::cout << "is-bez:" << std::boolalpha << geom::bspline_ops::is_bezier(p.first) << "\n";
  std::cout <<"at 0.3806:"<< p.second.eval(.3806) << "\n";
  std::cout <<"at 0.4:"<< p.second.eval(0.4) << "\n";
  std::cout <<"at 0.5:"<< p.second.eval(0.5) << "\n";
  std::cout <<"at 0.6:"<< p.second.eval(0.6) << "\n";
  std::cout <<"at 0.9:"<< p.second.eval(0.9) << "\n";

        std::cout << "---------------\nat 0:"<< p2.eval(0) << "\n";
  std::cout <<"at 0.2:"<< p2.eval(.2) << "\n";
  std::cout <<"at 0.3:"<< p2.eval(.3) << "\n";
  std::cout <<"at 0.3806:"<< p2.eval(.3806) << "\n";
  std::cout <<"at 0.4:"<< p2.eval(0.4) << "\n";
  std::cout <<"at 0.5:"<< p2.eval(0.5) << "\n";
  std::cout <<"at 0.6:"<< p2.eval(0.6) << "\n";
  std::cout <<"at 0.9:"<< p2.eval(0.9) << "\n";

   std::cout << "---------------\nat 0:"<< bs.eval(0) << "\n";
  std::cout <<"at 0.2:"<< bs.eval(.2) << "\n";
  std::cout <<"at 0.3:"<< bs.eval(.3) << "\n";
  std::cout <<"at 0.3806:"<< bs.eval(.3806) << "\n";
  std::cout <<"at 0.4:"<< bs.eval(0.4) << "\n";
  std::cout <<"at 0.5:"<< bs.eval(0.5) << "\n";
  std::cout <<"at 0.6:"<< bs.eval(0.6) << "\n";
  std::cout <<"at 0.9:"<< bs.eval(0.9) << "\n";
   }

{   cpts2d_t pts(4);
  pts[0] = make_pt(0.0,1);
  pts[1] = make_pt(0.4,1);
  pts[2] = make_pt(0.2,1);
  pts[3] = make_pt(-0.2, 1);
  //pts[4] = make_pt(0.0,0.0);
  std::vector<double> ks(5);
  ks[0] = 0.0; ks[1] = 0.3; ks[2] = 0.6;
  ks[3] = 0.8; ks[4] = 1.0;
  auto rbs = geom::make_rbspline(
      geom::make_periodic_bspline_wrap(pts,ks,2));
  std::cout << rbs.eval(0.1) << "\n";
}

{
    
      cpts3d_t pts(4);
  pts[0] = make_pt(0.0,0.0,1);
  pts[1] = make_pt(0.4,0.3,1);
  pts[2] = make_pt(0.2,0.8,1);
  pts[3] = make_pt(-0.2, 0.4,1);
  //pts[4] = make_pt(0.0,0.0);
  std::vector<double> ks(5);
  ks[0] = 0.0; ks[1] = 0.3; ks[2] = 0.6;
  ks[3] = 0.8; ks[4] = 1.0;
  geom::periodic_bspline<geom::point3d_t> bs(geom::make_periodic_bspline_wrap(pts, ks,2));
  geom::rational_bspline<geom::periodic_bspline<geom::point3d_t> >
      rbs = geom::make_rbspline(std::move(bs));
  auto const & pc = geom::bspline_ops::rotate_base_knot(rbs.spline(), 3);
  auto rpc = geom::bspline_ops::reparametrize(pc,0,1);
  
  auto rbs1 = geom::bspline_ops::split_periodic_curve(rbs, 0.1);
  rbs1.eval(0.1);
  auto const &vs = rbs1.eval_derivatives(2,0.1);
  std::cout << (vs.front()) << vs.back() << vs[1];

}
{
    geom::circle<geom::point2d_t> c( geom::make_pt(0.0,0.0),
        geom::make_pt(0,1.0) ); 
    std::cout << c.eval(0) << "," << c.eval(M_PI/2) <<"\n" << c.eval(M_PI)  << ", "<<  c.eval(3*M_PI/2) << "\n";
}
{
    geom::point2d_t pts[3];
    pts[0] = make_pt(0.0,0.0);
    pts[1] = make_pt(0.4,0.3);
    pts[2] = make_pt(1.0,0.0);
    geom::vector2d_t v[2];
    v[0] = geom::make_vec(0,1);
    v[1] = geom::make_vec(0.4,0.5);
    

    geom::conic_arc<geom::point2d_t> c = geom::make_conic_arc( pts, v) ;
    std::cout << c.eval(0) << "," << c.eval(0.5) << c.eval(1.0)<< "\n" ;
    std::cout << "type:" <<( (c.type()== geom::hyperbola )? "hyperbola\n"
        :(c.type()== geom::ellipse )? "ellipse\n"
        :(c.type()== geom::parabola )? "parabola\n" : "unknown\n") ;
    auto rbs ( geom::make_rbspline0(c));
    std::cout << rbs.eval(0) << "," << rbs.eval(0.5) << rbs.eval(1.0)<< "\n" ;
    
}
   
  return 0;
}

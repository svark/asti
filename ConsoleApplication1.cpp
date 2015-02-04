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

geom::point2d_t make_pt(double s,double t)
{
  double p[] = {s,t};
  return geom::point2d_t(p);
}

std::ostream& operator<<( std::ostream& os, const geom::point2d_t& pt)
{
  os << "(" <<pt.p[0] <<"," << pt.p[1] << ")";
  return os;
}

#include "geom_exception.hpp"
#include "bspline.hpp"
#include "bspline_ops.hpp"
double cube(double u) { double x = u ; return x*x*x -3*x+1;}

int _tmain(int argc, _TCHAR* argv[])
{

  std::vector<geom::point2d_t> pts(5);
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
  std::cout <<"diff at .20, 0.19:" <<  f0 + 0.01 * f1  + (0.01 * 0.01* f2)/2.0  + (0.01 * 0.01  *0.01 * f3)/6.0 - std::get<0>(fn(0.20)) << "\n";


  std::cout <<"diff at .21, 0.19" <<  geom::make_pt((bs.eval(.19) + 0.02 * bs.eval_derivative(1,.19) + 0.02 *0.01 *
                                                     bs.eval_derivative(2,0.19) )- bs.eval(0.21)) << "\n" ;
  std::cout <<"dot:" << dot(geom::make_vec( make_pt(1.0 , 0)) , geom::make_vec(make_pt(0.0,1.0) )) << "\n";

  std::cout << '(' << f0 <<',' <<  f1  << ')' << '\n';
  const double ts[] = {-1,-1,-1,-1,0.5,1,1,1,1};
  auto spl0 = geom::bspline_ops::cubic_approx1d(std::function<double(double)>(cube),std::vector<double>( ts , ts + sizeof(ts)/sizeof(ts[0]) ));
  std::cout << "cube approx at -1,0.1,0,1:" << spl0.eval(-1) << "," << spl0.eval(0.1) << ',' << spl0.eval(0) << ',' << spl0.eval(1) << "\n";

  std::vector<double> newknots(bs.knots().size() + 2);
  double arr[] = {0.225,0.75};
  std::merge(bs.knots().cbegin(), bs.knots().cend(), arr, arr + 2, newknots.begin());
  newknots.insert(newknots.begin(), 0.0);
  newknots.insert(newknots.end(), 1.0);
  auto spl2 = geom::bspline_ops::cubic_approx1d(
                                                std::function<double(double)>(v_vdash),
                                                newknots);
  try {
    auto pt_ = spl2.eval(-1);
    std::cout << "should not reach here\n";
  }catch(geom::geom_exception e) {
    std::cout << e.what() << "\n";
  }
  std:: cout << "\n dist.dist'  at 0.1,0.2,1:" << v_vdash(0.25) << ',' << v_vdash(0.2) << ',' << v_vdash(.5) << "\n";
  std:: cout << "\n dist.dist' approx at 0.1,0.2,1:" << spl2.eval(0.25) << ',' << spl2.eval(0.2) << ',' << spl2.eval(.5) << "\n";
  std::cout << "soln NEAR param .2::" << geom::bspline_ops::foot_param(bs,p);
  auto spls = geom::bspline_ops::split_into_bezier_patches(bs);
  
  double us[] = {0.1,0.6,0.9};
  int i=0;
  for(auto s : spls)
  {
	  std::cout << "bez:" << s.eval(us[i]) << "spl:" << bs.eval(us[i]) << "\n";
	  ++i;
  }
  //auto spls2 = geom::bspline_ops::split_into_bezier_patches_hard(bs);
  //bool ok = std::equal(spls.begin(), spls.end(),spls2.begin() );
  //std::cout << std::boolalpha << ok;
  return 0;
}

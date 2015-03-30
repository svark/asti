#ifndef CIRCLE_HPP
#define CIRCLE_HPP
#include "rational_bspline_cons.hpp"
#define TY(EXP) typename EXP::type
/*

*/
#include <algorithm>
namespace geom{
  template <class Point>
  struct circle
  {
      enum {dim = point_dim<Point>::dimension};
      typedef Point point_t;
      typedef decltype(make_vec<point_t>(point_t())) vector_t;

      template<class Point>
      struct allow_ydir_if_dim_3
      {
          typename std::enable_if< point_dim<PointU>::dimension == 3
                                   && dim ==3, vector3d_t const &>::type type;
      };

      template<class Point>
      struct allow_ydir_if_dim_3
      {
          typename std::enable_if< point_dim<PointU>::dimension == 3
                                   && dim ==3, vector3d_t const &>::type type;
      };

      template <class PointU>
      circle(const PointU& center_,
             const PointU& point_,
             allow_ydir_if_dim_3<PointU>::type ydir_)
          :center(center_), start_pt(point_), ydir(ydir_)
      {
          auto x = start_pt - center;
          assert(!tol::eq(len(x),0));
        ydir -= x*dot((start_pt-center),ydir)/len(x);
        ydir = normalize(ydir);
    }

     template <class PointU>
     circle(const PointU& center_,
           const PointU& point_,
           typename std::enable_if< point_dim<PointU>::dimension == 2
            && dim == 2, int >::type = 0 )
        :center(center_), start_pt(point_)
    {
        auto x = start_pt - center;
        ydir = normalize( make_vec(x[1],-x[0]) );
    }

    double
    static foot_param(const circle &c,
                      const point_t& p)
    {
      auto x = c.start_pt - c.center;
      x.normalize();
      auto y = ydir;
      y.normalize();

      auto v = (p - c.center);

      if( sqlen(v) < eps*eps)
        throw geom_exception(point_at_axis_error);

      return atan2(dot(v,y) *y , dot(v,x) * x);
    }

    template <class PointIter,class ParamIter>
    void
    static foot_param(const circle &c,
                      PointIter ps,
                      PointIter end, ParamIter out)
    {
      auto x = c.start_pt - c.center;
      x.normalize();
      auto y = ydir;
      y.normalize();

      for( ;ps!=end; ++ps, ++out ) {

        auto p =*ps;
        auto v = (p - c.center);

        if( sqlen(v) < eps*eps)
          throw geom_exception(point_at_axis_error);

        *out =  atan2(dot(v,y) *y , dot(v,x) * x);
      }
    }

    static circle<point_t>
    make_circle(const point_t& p1,
                const point_t& p2,
                const point_t& p3)
    {
      auto v1 = p1 - p2;
      const double eps = 1e-9;

      double lv1 = sqlen(v1);
      if(lv1 < eps*eps)
        throw geom_exception(circle_too_small);

      auto v2 = p2 - p3;
      double lv2 = sqlen(v2);
      if(lv2 < eps*eps)
        throw geom_exception(circle_too_small);

      auto v3 = p3 - p1;
      double lv3 = sqlen(v3);

      double d1 = sqlen( cross( v1, v2) );
      center = dot(v1 , -v3) / d*d;
      double lens[] ={lv1,lv2,lv3};
      auto vs[] = {v1,v2,v3};

      uint idx = std::max_element(lens,lens+3) - lens;

      v1 = vs[(idx + 1)%3];   lv1 = lens[(idx + 1)%3];
      v2 = vs[(idx + 2)%3];   lv2 = lens[(idx + 2)%3];
      v3 = vs[idx];           lv3 = lens[idx];

      point_t pts[] = {p1,p2,p3};
      std::rotate(pts, pts + (idx + 1)%3, pts + 3);
      auto normal = cross(v1,v2);
      double denom = sqlen(normal);
      if(denom < eps*eps )
        throw geom_exception(circle_too_small);

      double alpha = lv2 * dot(v1,-v3) *0.5 / ( denom );
      double beta  = lv3 * dot(-v1,v2) *0.5 / ( denom );

      return circle<point_t>(lerp(alpha,beta,pts[0],pts[1],pts[2]),
                         pts[0], cv/(-sqrt(denom)) );
    }
    //  (@file :file-name "media/circle2.png" :to "./media/circle2.png" :display "circle to nurbs")
    point_t eval(double u) const
    {
      au
      auto y = ydir;
      return center + x * cos(u) + y * sin(u);
    }

    template <class ParamIter ,class PointIter>
    void eval(ParamIter us, ParamIter end, PointIter out) const
    {
      auto x = (start_pt - center);
      auto y = ydir;
      for( ;us!=end;++us,++out)
         *out = center + x * cos(u) + y * sin(u);
    }

    vector_t tangent(double u) const
    {
      auto x = (start_pt - center);
      auto y = ydir;
      return  - x * sin(u) + y * cos(u);
    }

   vector_t normal(double u) const
    {
      return eval(u) - center;
    }

  private:
    point_t  center;
    point_t  start_pt;
    vector_t ydir;
  };

  template <class Point>
static rational_bspline< periodic_bspline<Point> >
    to_rational(const circle<Point>& c)
    {
      auto x = (c.start_pt - c.center);
      auto y = ydir;

      auto a = c.start_pt + y * 2 * c.radius  ;
      auto c = c.start_pt - y * c.radius +  x * 2 * cos(M_PI/3) * c.radius;
      auto b = c.center -  y * c.radius  - x * 2 * cos(M_PI/3) * c.radius;
      auto p = lerp(0.5,a,c);
      auto q = lerp(0.5,a,b);
      auto r = lerp(0.5,b,c);

      double weights[] = {1,0.5,1,0.5,1,0.5}
      point_t cpts[] = {p,  a   ,q    ,b    ,r        ,c/*, p*/};
      double ts[]      = {0,1/6, 1/3,1/3, 0.5, 2/3, 2/3,5/6, 1};
      return make_periodic_rbspline(cpts,weights,ts,2);
    }

}




#endif

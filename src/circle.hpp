#ifndef CIRCLE_HPP
#define CIRCLE_HPP
/*
(@file :file-name "circle.png" :to "./circle.png" :display "circle")
*/
#include <algorithm>
namespace geom{
  template <dim>
  struct circle
  {
    circle(const pt_t<dim>& center_,
           const pt_t<dim>& point_,
           const vec_t<dim> &normal_)
      :center(center_), start_point(point_), normal(normal_)
    {

    }

    double
    static foot_param(const circle &c,
                      const pt_t<dim>& p)
    {
      auto v = (p - c.center);
      v -= dot(v,c.normal)*c.normal;
      if( sqlen(v) < eps*eps)
        throw geom_exception(point_at_axis_error);

      auto x = c.start_pt - c.center;
      x.normalize();
      auto y = cross(c.normal,x);
      y.normalize();
      return atan2(dot(v,y) *y , dot(v,x) * x);
    }

    void
    static foot_param(const circle &c,
                      PointIter ps, PointIter end, ParamIter out)
    {
      auto x = c.start_pt - c.center;
      x.normalize();
      auto y = cross(c.normal,x);
      y.normalize();

      for( ;ps!=end; ++ps, ++out ) {

        auto p =*ps;
        auto v = (p - c.center);
        v -= dot(v,c.normal)*c.normal;

        if( sqlen(v) < eps*eps)
          throw geom_exception(point_at_axis_error);

        *out =  atan2(dot(v,y) *y , dot(v,x) * x);
      }
    }

    static circle<dim>
    make_circle(const pt_t<dim>& p1,
                const pt_t<dim>& p2,
                const pt_t<dim>& p3)
    {
      vec_t<dim> v1 = p1 - p2;
      const double eps = 1e-9;

      double lv1 = sqlen(v1);
      if(lv1 < eps*eps)
        throw geom_exception(circle_too_small);

      vec_t<dim> v2 = p2 - p3;
      double lv2 = sqlen(v2);
      if(lv2 < eps*eps)
        throw geom_exception(circle_too_small);

      vec_t<dim> v3 = p3 - p1;
      double lv3 = sqlen(v3);

      double d1 = sqlen( cross( v1, v2) );
      center = dot(v1 , -v3) / d*d;
      double lens[] ={lv1,lv2,lv3};
      vec_t<dim> vs = {v1,v2,v3};

      uint idx = std::max_element(lens,lens+3) - lens;
      v1 = vs[(idx + 1)%3]; lv1 = lens[(idx + 1)%3];
      v2 = vs[(idx + 2)%3]; lv2 = lens[(idx + 2)%3];
      v3 = vs[idx]; lv3 = lens[idx];
      pt_t<dim> pts[] = {p1,p2,p3};
      std::rotate(pts, pts + (idx + 1)%3, pts + 3);
      vec_t<dim> normal = cross(v1,v2);
      double denom = sqlen(normal);
      if(denom < eps*eps )
        throw geom_exception(circle_too_small);

      double alpha = lv2 * dot(v1,-v3) *0.5 / ( denom );
      double beta  = lv3 * dot(-v1,v2) *0.5 / ( denom );

      return circle<dim>(lerp(alpha,beta,pts[0],pts[1],pts[2]),
                         pts[0], cv/(-sqrt(denom)) );
    }
(@file :file-name "circle2.png" :to "./circle2.png" :display "circle to nurbs")
    static rational_bspline<dim>
    to_rational(const circle<dim>& c)
    {
      vec_t<dim> z = c.normal;
      vec_t<dim> x = (c.start_pt - c.center);
      vec_t<dim> y = cross(z,x);

      auto a = c.start_pt + c.radius/tan(M_PI/6) * y;
      auto c = c.start_pt - c.radius/tan(M_PI/6) * y;
      auto b = c.center - 2 * x;
      auto p = lerp(0.5,a,c);
      auto q = lerp(0.5,a,b);
      auto r = lerp(0.5,b,c);

      double weights[] = {1,0.5,1,0.5,1,0.5}
      pt_t<dim> cpts[] = {p,a,q,b,r,c,p};
      double ts[]      = {0,0,0,1/3,1/3,2/3,2/3,1,1,1};
      return rational_bspline(cpts,weights,ts,2);
    }

    pt_t<dim> eval(double u) const
    {
      vec_t<dim> z = normal;
      vec_t<dim> x = (start_pt - center);
      vec_t<dim> y = cross(z,x);
      return center + x * cos(u) + y * sin(u);
    }

    template <class ParamIter >
    void eval(ParamIter us, ParamIter end, PointIter out) const
    {
      vec_t<dim> z = normal;
      vec_t<dim> x = (start_pt - center);
      vec_t<dim> y = cross(z,x);
      for( ;us!=end;++us,++out)
         *out = center + x * cos(u) + y * sin(u);
    }

    vec_t<dim> tangent(double u) const
    {
      vec_t<dim> z = normal;
      vec_t<dim> x = (start_pt - center);
      vec_t<dim> y = cross(z,x);
      return  - x * sin(u) + y * cos(u);
    }

   vec_t<dim> normal(double u) const
    {
      return eval(u) - center;
    }

  private:
    pt_t<dim>  center;
    pt_t<dim>  start_pt;
    vec_t<dim> normal;
  };
}

#endif

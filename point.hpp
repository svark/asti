#ifndef _POINT_HPP
#define _POINT_HPP
#include <Eigen/Core>
#include <vector>
namespace geom {

  template <int dim>
  struct pt_t {
  
	  pt_t(){}
      
      pt_t(const double v[]): p(v)
      {
      }

      explicit pt_t(double v)
      {
          for(int i = 0; i < dim;++i)
              p(i) = v;
      }
	  
	  pt_t(Eigen::Matrix<double,dim,1>&& v) {p.swap(v);}

	  pt_t(const Eigen::Matrix<double,dim,1>& v):p(v) {}

      double operator[](std::size_t i) const { return p[i];}
      double& operator[](std::size_t i) { return p[i];}
	   	 
	  template <class VecT>
	  pt_t& operator+=( const VecT& v) { p+=v.v; return *this;}
	  template <class VecT>
	  pt_t& operator-=( const VecT& v) { p-=v.v; return *this;}

	  pt_t& operator+=(const double& v) { for(int i =0; i < dim;++i) p[i]+=v; return *this;}
	  template <class VecT>
	  pt_t& operator-=(const double& v) { for(int i =0; i < dim;++i) p[i]-=v; return *this;}
	  
      Eigen::Matrix<double,dim,1>  p;
      enum { dimension = dim };
  };

  template struct pt_t<1>;
  template struct pt_t<2>;
  template struct pt_t<3>;
  template struct pt_t<4>;

  template <class Point>
  struct point_dim
  {
	  enum {dimension = Point::dimension};
  };
  
  template <>
  struct point_dim<double>
  {
	  enum {dimension = 1};
  };

  template <class PointIter>
  struct point_traits
  {
      typedef typename std::iterator_traits<PointIter>::value_type point_t;
	  enum {dim = point_dim<point_t>::dimension};
  };
  

  typedef pt_t<1> point1d_t;
  typedef pt_t<2> point2d_t;
  typedef pt_t<3> point3d_t;
  typedef pt_t<4> point4d_t;

  template <class PointIter>
  static typename point_traits<PointIter>::point_t 
  centroid (PointIter pts,PointIter end)
  {
    typedef typename point_traits<PointIter>::point_t point_t;

    point_t avg = point_t();
	point_t zero(0.0);
    int num_pts = 0;
    for(;pts!=end;++pts,++num_pts) {
        avg += (*pts - zero);
    }
    scale(avg , 1.0/num_pts);
    return avg;
  }

  template<int dim>
  static pt_t<dim> zero_pt()
  {
    pt_t<dim> pt(0);
    return pt;
  }
  template <int dim>
  bool operator!=(const pt_t<dim>& p1, const pt_t<dim>& p2) { return !(p1.p==p2.p); }
  template <int dim>
  bool operator==(const pt_t<dim>& p1, const pt_t<dim>& p2) { return sqlen(p1-p2)<1e-16; }

  template <int dim>
  struct vec_t
  {
	  vec_t() {}

      vec_t(const double ve[]): v(ve)
      {
      }
	  vec_t(Eigen::Matrix<double,dim,1>&& w):v(w)
	  {
	  }
	  vec_t(const Eigen::Matrix<double,dim,1>& w):v(w)
	  {
	  }

      explicit vec_t(double ve)
      {
          for(int i = 0; i < dim;++i)
              v[i] = ve;
      }

	  template <class VecT>
	  vec_t& operator+=( const VecT& v) { this->v+=v.v; return *this;}
	  template <class VecT>
	  vec_t& operator-=( const VecT& v) { this->v-=v.v; return *this;}

      vec_t& negate() { v = -v; return *this;}
      double operator[](std::size_t i) const { return v[i];}
      double& operator[](std::size_t i) { return v[i];}
      Eigen::Matrix<double,dim,1> v;
	  enum { dimension = dim };
  };
  
  template <int point_dim>
  bool operator==(const vec_t<point_dim>&v1,const vec_t<point_dim>&v2) { return sqlen(v1-v2) < 1e-16;}
  

  template struct vec_t<2>;
  template struct vec_t<3>;
  template struct vec_t<4>;

  template<int dim>
  vec_t<dim> zero_vec()
  {
    vec_t<dim> vec;
    for(int i = 0; i < dim; ++i)
      vec.v[i] = 0.0;
    return vec;
  }

  template<int dim>
  double len(const vec_t < dim >& vec)
  {
      return vec.v.norm();
  }

  template<int dim>
  double odist(pt_t< dim >& pt)
  {
      return pt.p.norm();
  }

  template<class point_t,int dim>
  point_t operator+(const point_t & pt, const vec_t<dim>& vec)
  {
    auto res = pt;
    static_assert(point_t::dimension == dim,
                  "dimensions of point and vector "
                  "must match for addition");
    res += vec;
    return res;
  }

  template<class point_t,int dim>
  point_t operator-(const point_t & pt,
                    const vec_t<dim>& vec)
  {
    auto res = pt;
    static_assert(point_t::dimension == dim,"dimensions of point and vector must match for subtraction");
    res -= vec;
    return res;
  }
  template <int dim>
  vec_t<dim> operator / ( vec_t<dim> &vec, double s)
  {
    return vec.v*(1/s);
  }

  template <int dim>
  vec_t<dim> operator * ( vec_t<dim> &vec, double s)
  {
    return vec.v * s;
  }
  template <int dim>
  vec_t<dim> operator / ( double s, vec_t<dim> &vec )
  {
    return vec.v*(1/s);
  }

  template <int dim>
  vec_t<dim> operator * ( double s, vec_t<dim> &vec)
  {
    return vec.v * s;
  }


  template <int dim>
  pt_t<dim> copy_pt(const point2d_t& p_)
  {
    pt_t<dim> pt;
    static_assert(dim>=2,"atleast 2 dim");
    pt.p[0] = p_.p[0];
    pt.p[1] = p_.p[1];
    for(int j = 2;j < dim;++j)
      pt.p[j] = 0;
    return pt;
  }

  template <int dim>
  pt_t<dim> copy_pt(const point3d_t& p_)
  {
    pt_t<dim> pt;
    static_assert(dim>=2,"atleast 2 dim");
    int low = std::min(dim,3);
    for(int j = 0;j < low;++j)
      pt.p[j] = p_.p[j];
    for(int j = low;j < dim;++j)
      pt.p[j] = 0.0;
    return pt;
  }

  template <int dim>
  pt_t<dim> copy_pt(const point4d_t& p_)
  {
    pt_t<dim> pt;
    static_assert(dim>=2,"atleast 2 dim");
    int low = std::min(dim,4);
    for(int j = 0;j < low;++j)
      pt.p[j] = p_.p[j];
    for(int j = low;j < dim;++j)
      pt.p[j] = 0.0;
    return pt;
  }

  template <int dim>
  pt_t<dim> make_pt(const vec_t<dim>& v)
  {
    pt_t<dim> pt;
    pt.p = v.v;
    return pt;
  }
  template <int dim>
  pt_t<dim> make_pt(vec_t<dim>&& v)
  {
	return pt_t<dim>(std::forward<Eigen::Matrix<double,dim,1> >(v.v));
  }

  template <class point_t>
  vec_t<point_dim<point_t>::dimension> make_vec(const point_t & pt)
  {
    return vec_t<point_dim<point_t>::dimension>(pt.p);
  }
  template <class point_t>
  vec_t<point_dim<point_t>::dimension> make_vec( point_t&& pt)
  {
    enum { dim = point_dim<point_t>::dimension };
    return vec_t<dim>(
		std::forward<Eigen::Matrix<double,dim,1> > (
			pt.p));
  }

  inline double make_vec(double  pt)
  {
    return pt;
  }

  typedef vec_t<1> vector1d_t;
  typedef vec_t<2> vector2d_t;
  typedef vec_t<3> vector3d_t;
  typedef vec_t<4> vector4d_t;

  template <int dim>
  vec_t<dim>
  operator-(const pt_t<dim>& p1, const pt_t<dim>& p2)
  {
    vec_t<dim> vec;
    vec.v = p1.p - p2.p;
    return vec;
  }

  //{{{  (@* "3d vector ops")
  inline vec_t<3>
  cross(const vec_t<3> &v1, const vec_t<3>& v2)
  {
    vec_t<3> v3;
    v3.v = v1.v.cross(v2.v); // todo: move constructors?
    return v3;
  }

  inline vec_t<4>
  cross(const vec_t<4> &v1, const vec_t<4>& v2)
  {
    vec_t<4> v3;
    v3.v[0] = (v1.v[1]*v2.v[2] - v1.v[2]*v2.v[1]);
    v3.v[1] = (v1.v[2]*v2.v[0] - v1.v[0]*v2.v[2]);
    v3.v[2] = (v1.v[0]*v2.v[1] - v1.v[1]*v2.v[0]);
    v3.v[3] = v1.v[3]*v2.v[3];
    return v3;
  }

  inline vec_t<3>
  cross(const vec_t<2> &v1, const vec_t<2>& v2)
  {
    vec_t<3> v3;
    v3.v[0] = 0;
    v3.v[1] = 0;
    v3.v[2] = (v1.v[0]*v2.v[1] - v1.v[1]*v2.v[0]);
    return v3;
  }

  template <int dim>
  double dot(const vec_t<dim>&v1,
             const vec_t<dim>&v2)
  {
    return v1.v.dot(v2.v);
  }
  
  template <int dim>
  double angle_between(const vec_t<dim>&v1,
                       const vec_t<dim>&v2)
  {
      return atan2(len(cross(v1,v2)),dot(v1,v2));
  }

// effectively a complicated looking macro
 template <class T>
 decltype(t[j+1] - t[j]) E(T t,std::size_t j, std::size_t n = 1)
 {
     return t[j+n] - t[j];
 }
 
 template <int dim>
 double angle_between(const vec_t<dim>&v1,
                      const vec_t<dim>&v2, 
                      const vec_t<dim>& axis)
  {
      auto nrml = cross(v1,v2);
      auto sign = dot(nrml,axis) > 0 ? 1 : -1;
      return atan2(len(nrml),dot(v1,v2))*sign;
  }

  template <int dim>
  double sqlen(const vec_t<dim>& v)
  {
     return v.v.squaredNorm();
  }

  //}}}
  //{{{ (@* "linear interpolation methods")
 
  extern point4d_t
  lerp( double lambda,
        const point4d_t& p1,
        const point4d_t& p2);
  
  extern point3d_t
  lerp( double lambda ,
        const point3d_t& p1,
        const point3d_t& p2);

  extern point2d_t
  lerp( double lambda ,
        const point2d_t& p1,
        const point2d_t& p2);

  extern double
  lerp( double lambda ,
        const double& p1,//0
        const double& p2);//1

  ////////////////////////////////
  extern point4d_t
  lerp( double lambda , double mu,
        const point4d_t& p1,//1,0
        const point4d_t& p2,//0,1
        const point4d_t& p3//0,0
        );

  extern point3d_t
  lerp( double lambda , double mu,
        const point3d_t& p1,
        const point3d_t& p2,
        const point3d_t& p3
        );
 
  extern point2d_t
  lerp( double lambda , double mu,
        const point2d_t& p1,
        const point2d_t& p2,
        const point2d_t& p3
        );

  template <int dim>
  pt_t<dim>
  dlerp( double lambda, 
		const pt_t<dim> &p1, const pt_t<dim>&p2) 
  { 
	return pt_t<dim>((p2.p - p1.p)*lambda);
  }
  
  inline double  
  dlerp( double lambda, 
		const double &p1, const double&p2) 
  { 
	return ((p2 - p1)*lambda);
  }

  template <int dim>
  void  scale(pt_t<dim>& p, double fac) {
      p.p*=fac;
  }

  inline double scale(double v, double fac) {
      return v*fac;
  }
  inline double sqlen(double v ) {
	  return v*v;
  }

  inline double dot(double v , double w) {
	  return v*w;
  }

  template <int dim>
  pt_t<dim>  scaled_copy(pt_t<dim> p, double fac) {
      p.p*=fac;
      return p;
  }
  template <int dim>
  const vec_t<dim>& operator - ( vec_t<dim>&& v ) 
  {
	  vec_t<dim> w(-v.v);
	  return v;
  }

  template <int dim>
  const Eigen::Matrix<double,dim,1>& eigen_vec(const vec_t<dim>& ve)
  {
      return ve.v;
  }
  //}}}

  //{{{  interleave points with weights
  template<int dim>
  std::vector<pt_t<dim+1>>
  interleave(const std::vector<pt_t<dim>> &vs,
             const std::vector<double>  &ws)
  {
    std::vector<pt_t<dim+1> > vws;
    vws.reserve(vs.size());
    auto witer = ws.begin();
    for(const pt_t<dim> &v: vs)
      {
        pt_t<dim+1> vw = ( copy_pt<dim+1>(v) );
        vw.p[dim] = *witer;
        ++witer;
        vws.push_back(vw);
      }
    return vws;
  }
  //}}}
}
namespace std
{		
	template<int dim>
	void swap( geom::pt_t<dim>& p1, geom::pt_t<dim>&p2)
    {
		p1.p.swap(p2.p);
    }
	template<int dim>
	void swap( geom::vec_t<dim>& v1, geom::vec_t<dim>&v2)
    {
		v1.v.swap(v2.v);
    }
    inline double fmin(double a, double b ) { return a < b? a: b;}
    inline double fmax(double a, double b ) { return b < a? a: b;}
}
#endif//_POINT_HPP

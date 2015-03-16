#ifndef _POINT_HPP
#define _POINT_HPP
#include <Eigen/Core>
#include "tol.hpp"
#include <type_traits>
namespace geom {

template <int dim>
struct pt_t {

    typedef Eigen::Matrix<double,dim,1> EMT;

    pt_t(){}

    pt_t(const double v[]): p(v)
    {
    }

    pt_t(const pt_t& other): p(other.p)
    {
    }

    template <int odim>
    pt_t(const pt_t<odim>& other, double fill = 0.0)
    {
        const int low = std::min(dim,odim);
        for(int j = 0;j < low;++j)
            p[j] = other.p[j];
        for(int j = low;j < dim;++j)
            p[j] = fill;
    }


    explicit pt_t(double v)
    {
        for(int i = 0; i < dim;++i)
            p(i) = v;
    }

    pt_t(const EMT& v)
        :p(v)
    {}


    template <class Derived,int BlockRows, int BlockCols, bool InnerPanel>
    pt_t(const Eigen::Block<Derived,BlockRows,BlockCols,InnerPanel>& b)
        :p(EMT(b))
    {
    }

    double operator[](std::size_t i) const { return p[i];}
    double& operator[](std::size_t i) { return p[i];}

    template <class VecT>
    pt_t& operator+=( const VecT& v) { p+=v.v; return *this;}
    template <class VecT>
    pt_t& operator-=( const VecT& v) { p-=v.v; return *this;}

    pt_t& operator+=(const double& v) {
        for(int i =0; i < dim;++i) p[i]+=v; return *this;}

    pt_t& operator-=(const double& v) {
        for(int i =0; i < dim;++i) p[i]-=v; return *this;}


    pt_t& operator*=(const double& v) {
        for(int i =0; i < dim;++i) p[i]*=v; return *this;}

    pt_t& operator/=(const double& v) {
        for(int i =0; i < dim;++i) p[i]/=v; return *this;}

    EMT& get() { return p; }
    const EMT& cget() const { return p; }

    EMT  p;
    enum { dimension = dim };
    enum { NeedsToAlign = (sizeof(EMT)%16)==0 };
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign)
};

template struct pt_t<1>;
template struct pt_t<2>;
template struct pt_t<3>;
template struct pt_t<4>;


typedef pt_t<1> point1d_t;
typedef EIGEN_ALIGN16 pt_t<2> point2d_t;
typedef pt_t<3> point3d_t;
typedef EIGEN_ALIGN16 pt_t<4> point4d_t;

template <class PointIter>
auto
centroid (PointIter pts,PointIter end) ->  typename std::decay<decltype(pts[0])>::type
{
    typedef typename std::decay<decltype(pts[0])>::type point_t;

    point_t avg = point_t();
    point_t zero(0.0);
    int num_pts = 0;
    for(;pts!=end;++pts,++num_pts) {
        avg += (*pts - zero);
    }
    scale(avg , 1.0/num_pts);
    return avg;
}

template <int dim>
bool operator!=(const pt_t<dim>& p1, const pt_t<dim>& p2) {
    return !(p1.p==p2.p);
}
template <int dim>
bool operator==(const pt_t<dim>& p1, const pt_t<dim>& p2) {
    return sqlen(p1-p2)<1e-16;
}


template <int dim>
struct vec_t
{
    enum { NeedsToAlign = (sizeof(Eigen::Matrix<double,dim,1>)%16)==0 };
    typedef Eigen::Matrix<double,dim,1> EMT;

    vec_t() {}

    vec_t(const double ve[]): v(ve)
    {
    }
    vec_t(const vec_t&other):v(other.v){}

    vec_t(const EMT& w):v(w)
    {
    }

    template <int odim>
    vec_t(const vec_t<odim>& other)
    {
        const int low = std::min(dim,2);
        for(int j = 0;j < low;++j)
            v[j] = other[j];
        for(int j = low;j < dim;++j)
            v[j] = 0.0;
    }

    explicit vec_t(double ve) {
        for(int i = 0; i < dim;++i)
            v[i] = ve;
    }

    vec_t& operator+=( const vec_t& v) { this->v+=v.v; return *this;}
    vec_t& operator-=( const vec_t& v) { this->v-=v.v; return *this;}

    vec_t& operator*=( double a) { this->v*=a; return *this;}
    vec_t& operator/=( double a) { this->v/=a; return *this;}

    vec_t& negate() { v = -v; return *this;}
    double operator[](std::size_t i) const { return v[i];}
    double& operator[](std::size_t i) { return v[i];}
    enum { dimension = dim };

    EMT& get() { return v; }
    const EMT& cget() const { return v; }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign)
    EMT v;
};

template <int dim>
bool operator==(const vec_t<dim>&v1,const vec_t<dim>&v2) {
    return sqlen(v1-v2) < 1e-16;
}




template struct vec_t<2>;
template struct vec_t<3>;
template struct vec_t<4>;


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
    point_t res(pt);
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
    point_t res(pt);
    static_assert(point_t::dimension == dim,
                  "dimensions of point and vector "
                  "must match for subtraction");
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


inline point3d_t  lower_dim(const point4d_t& pt)
{
    return point3d_t(pt);
}

inline point2d_t  lower_dim(const point3d_t& pt)
{
    return point2d_t(pt);
}

inline double lower_dim(const point2d_t& pt)
{
    return pt[0];
}




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
auto E(T t,std::size_t j, std::size_t n = 1)->decltype(t[j+1] - t[j])
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

template <class XprType,int BlockRows,int BlockCols,bool InnerPanel>
Eigen::Block<const XprType,BlockRows,BlockCols,InnerPanel>
lerp( double lambda,
      const Eigen::Block<XprType,BlockRows,BlockCols,InnerPanel>& p1,
      const Eigen::Block<XprType,BlockRows,BlockCols,InnerPanel>& p2)
{
   return Eigen::Block<const XprType,BlockRows,BlockCols,InnerPanel>(
       ((1-lambda)*p1+lambda*p2).eval(),0,0);
}
/*
Eigen::Block<Eigen::Matrix<double,-1,1>,1,1,false>
lerp( double lambda,
      const Eigen::Block<Eigen::Matrix<double,-1,1>,1,1,false>& p1,
      const Eigen::Block<Eigen::Matrix<double,-1,1>,1,1,false>& p2)
{
    Eigen::Block<Eigen::Matrix<double,-1,1>,1,1,false> tmp(p1);
    tmp[0] = lerp(lambda,p1[0],p2[0]) ;
    return tmp;
}
*/
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
    return (p2 - p1)*lambda;
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

inline double len(double v ) {
    return fabs(v);
}

inline double dot(double v , double w) {
    return v*w;
}

template <class point_t>
point_t  scaled_copy( const point_t&  p, double fac) {
    point_t q(p);
    q*=fac;
    return q;
}

template <int dim>
const vec_t<dim>& operator - ( vec_t<dim>&& v )
{
    vec_t<dim> w(-v.v);
    return v;
}

template <class PointVec>
auto eigen_vec(const PointVec& ve) -> decltype(ve.cget())
{
    return ve.cget();
}

inline Eigen::Matrix<double,1,1> eigen_vec(double ve)
{
    Eigen::Matrix<double,1,1> v;
    v[0]=ve;
    return v;
}
//}}}



inline point2d_t make_pt(double s,double t)
{
  double p[] = {s,t};
  return point2d_t(p);
}


inline point3d_t make_pt(double s,double t, double w)
{
  double p[] = {s,t,w};
  return point3d_t(p);
}

inline point4d_t make_pt(double s,double t, double w, double x)
{
  double p[] = {s,t,w,x};
  return point4d_t(p);
}

typedef vec_t<1> vector1d_t;
typedef EIGEN_ALIGN16 vec_t<2> vector2d_t;
typedef vec_t<3> vector3d_t;
typedef EIGEN_ALIGN16 vec_t<4> vector4d_t;

inline vector3d_t lower_dim(const vector4d_t& v)
{
    return vector3d_t(v);
}

inline vector2d_t lower_dim(const vector3d_t& v)
{
    return vector2d_t(v);
}

inline double lower_dim(const vector2d_t& v)
{
    return v[0];
}

inline vector2d_t normalize(const vector2d_t& vec)
{
    vector2d_t v( vec );
    v *= (1.0/v.v.norm());
    return v;
}
inline vector3d_t normalize(const vector3d_t& vec)
{
    vector3d_t v( vec );
    v *= (1.0/v.v.norm());
    return v;
}

inline  vector4d_t normalize(const vector4d_t& vec)
{
    vector4d_t v( vec );
    v *= (1.0/v.v.norm());
    return v;
}


inline double normalize(const double & vec)
{
    return 1;
}




//{{{  interleave points with weights


inline vector2d_t _make_vec(const point2d_t & pt)
{
    return vector2d_t(pt.cget());
}

inline vector3d_t _make_vec(const point3d_t & pt)
{
    return vector3d_t(pt.cget());
}

inline vector4d_t _make_vec(const point4d_t & pt)
{
    return vector4d_t(pt.cget());
}

inline double _make_vec(const double& pt)
{
    return pt;
}

template <class Point>
auto make_vec( const Point &p ) -> decltype(_make_vec(p))  { return  _make_vec(p); } 


inline point2d_t _make_pt(const  vector2d_t& pt)
{
    return point2d_t(pt.cget());
}

inline point3d_t _make_pt(const  vector3d_t & pt)
{
    return point3d_t(pt.cget());
}

inline point4d_t _make_pt(const   vector4d_t& pt)
{
    return point4d_t(pt.cget());
}

inline double _make_pt(const double& v)
{
    return v;
}

template <class Vector>
auto make_pt( const Vector &p ) -> decltype(_make_pt(p))  { return  _make_pt(p); } 




inline vector2d_t make_vec(double s,double t)
{
  double p[] = {s,t};
  return vector2d_t(p);
}


inline vector3d_t make_vec(double s,double t, double w)
{
  double p[] = {s,t,w};
  return vector3d_t(p);
}

inline vector4d_t make_vec(double s,double t, double w, double x)
{
  double p[] = {s,t,w,x};
  return vector4d_t(p);
}


template <class Point>
struct point_dim
{
    
    enum {dimension = Point::dimension};
    typedef Eigen::aligned_allocator<Point> alloc_t;
  
};

template <>
struct point_dim<double>
{
    enum {dimension = 1};
    typedef std::allocator<double> alloc_t;
};


template <class point_t, int dim>
struct inc_dimension_helper
{
    enum {NeedsToAlign = pt_t < dim + 1 >::NeedsToAlign };
    
    typedef typename std::conditional<NeedsToAlign, EIGEN_ALIGN16 pt_t < dim + 1 >, 
        pt_t < dim + 1 > >::type    type;
};

template <int dim>
struct inc_dimension_helper < vec_t < dim >, dim >
{
    enum {NeedsToAlign = vec_t < dim + 1 >::NeedsToAlign };
    typedef typename std::conditional<NeedsToAlign, EIGEN_ALIGN16 pt_t < dim + 1 >, 
        vec_t < dim + 1 > >::type    type;
};

template<class PointVec>
struct inc_dimension
{
    enum { dimension =  point_dim < PointVec >::dimension };
    typedef typename inc_dimension_helper < PointVec, dimension
                                          >::type type;
};


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

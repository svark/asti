#ifndef ASTI_POINT_FWD_HPP
#define ASTI_POINT_FWD_HPP
#include <Eigen/Core>
#include "type_utils.hpp"
#include <vector>
namespace geom {
//{{{  (@* "Point and vector decls")
template <int dim> struct pt_t;
template <int dim> struct vec_t;

typedef pt_t<1> point1d_t;
typedef EIGEN_ALIGN16 pt_t<2> point2d_t;
typedef pt_t<3> point3d_t;
typedef EIGEN_ALIGN16 pt_t<4> point4d_t;

typedef EIGEN_ALIGN16 vec_t<2> vector2d_t;
typedef vec_t<3> vector3d_t;
typedef EIGEN_ALIGN16 vec_t<4> vector4d_t;
//}}}
//{{{  (@* "lower / higher dim points / vectors")
vector3d_t lower_dim(const vector4d_t& v);
vector2d_t lower_dim(const vector3d_t& v);
double     lower_dim(const vector2d_t& v);
vector4d_t higher_dim(const vector3d_t& v);
vector3d_t higher_dim(const vector2d_t& v);

point3d_t  lower_dim(const point4d_t& pt);
point2d_t  lower_dim(const point3d_t& pt);
double     lower_dim(const point2d_t& pt);

point4d_t  higher_dim(const point3d_t& pt);
point3d_t  higher_dim(const point2d_t& pt);
point2d_t higher_dim(double pt);
//}}}
//{{{  (@* "mk_stdvec from a single/array of points or vectors")

template <int dim>
STDVEC(vec_t<dim>)
    mk_stdvec(const vec_t<dim> &v);

template <int dim>
STDVEC(pt_t<dim>)
    mk_stdvec(const pt_t<dim> &v);

std::vector<double>
mk_stdvec(const double &v);

template <int dim>
STDVEC(vec_t<dim>)
    mk_stdvec(const vec_t<dim> * vb, const vec_t<dim> * ve);

template <int dim>
STDVEC(pt_t<dim>)
    mk_stdvec(const pt_t<dim> * vb, const pt_t<dim> * ve);

std::vector<double>
mk_stdvec(const double *vb, const double *ve);
//}}}
//{{{  (@* "linear interpolation methods")
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

template <class PointVec>
PointVec dlerp( double lambda , const PointVec&p1, const PointVec&p2);
double dlerp(double lambda, const double &, const double & );
//}}}
//{{{  (@* "convert vector to point and vice versa")
vector2d_t _make_vec(const point2d_t & pt);
vector3d_t _make_vec(const point3d_t & pt);
vector4d_t _make_vec(const point4d_t & pt);
double _make_vec(const double& pt);

template <class Point>
auto make_vec( const Point &p ) -> decltype(_make_vec(p));

point2d_t _make_pt(const  vector2d_t& pt);
point3d_t _make_pt(const  vector3d_t & pt);
point4d_t _make_pt(const  vector4d_t & pt);
double _make_pt(const double& pt);

template <class Point>
auto make_pt( const Point &p ) -> decltype(_make_pt(p));
//}}}
//{{{  (@* "dot and len methods on vectors")
template <int dim>
double dot(const vec_t<dim>&v1, const vec_t<dim>&v2);
double dot(double v , double w);

template <int dim>
double sqlen(const vec_t<dim>& v);
double sqlen(double v ) ;

template <int dim>
double len(const vec_t<dim>&v);
double len(double v );


//}}}
//{{{  (@* "convert to eigen vec")
Eigen::Matrix<double,1,1> eigen_vec(double ve);

template <int dim>
Eigen::Matrix<double, dim , 1>
eigen_vec(const pt_t<dim>& ve);

template <int dim>
Eigen::Matrix<double, dim , 1>
eigen_vec(const vec_t<dim>& ve) ;
//}}}
//{{{  (@* "misc methods")
template <class PointIter>
auto
centroid (PointIter pts,PointIter end) -> RAWTYPE(pts[0]);
//}}}

//{{{
extern std::integral_constant<int,0> X;
extern std::integral_constant<int,1> Y;
extern std::integral_constant<int,2> Z;
extern std::integral_constant<int,3> W;

template <class Point,int i>
double coord(const Point&, std::integral_constant<int, i>);

template <class PtVec>
void set_quiet_NaN(PtVec& pv);

inline void set_quiet_NaN(double &d);

template <int dim>
double plen(const pt_t<dim>& p);
inline double plen(double p );


//}}}
}

#define VECTOR_TYPE(PointU) RAWTYPE(geom::make_vec(PointU()))
#define ARRAY_TYPE(PointU) RAWTYPE(geom::mk_stdvec(PointU()))
#endif

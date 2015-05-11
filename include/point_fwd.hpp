#include <Eigen/Core>
#include <vector>
namespace geom {

template <int dim> struct pt_t;
template <int dim> struct vec_t;

typedef pt_t<1> point1d_t;
typedef EIGEN_ALIGN16 pt_t<2> point2d_t;
typedef pt_t<3> point3d_t;
typedef EIGEN_ALIGN16 pt_t<4> point4d_t;

typedef EIGEN_ALIGN16 vec_t<2> vector2d_t;
typedef vec_t<3> vector3d_t;
typedef EIGEN_ALIGN16 vec_t<4> vector4d_t;



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

vector2d_t _make_vec(const point2d_t & pt);
vector3d_t _make_vec(const point3d_t & pt);
vector4d_t _make_vec(const point4d_t & pt);
double _make_vec(const double& pt);

template <class Point>
auto make_vec( const Point &p ) -> decltype(_make_vec(p));

template <class Point>
auto make_vec( const Point &p )->decltype(_make_vec(p));

std::vector<vector2d_t,Eigen::aligned_allocator<vector2d_t>>
    mk_stdvec(const vector2d_t &v);
std::vector<vector3d_t,Eigen::aligned_allocator<vector3d_t>>
    mk_stdvec(const vector3d_t &v);
std::vector<vector4d_t,Eigen::aligned_allocator<vector4d_t>>
    mk_stdvec(const vector4d_t &v);
std::vector<point2d_t,Eigen::aligned_allocator<point2d_t>>
    mk_stdvec(const point2d_t &v);
std::vector<point3d_t,Eigen::aligned_allocator<point3d_t>>
    mk_stdvec(const point3d_t &v);
std::vector<point4d_t,Eigen::aligned_allocator<point4d_t>>
    mk_stdvec(const point4d_t &v);
std::vector<double>
mk_stdvec(const double &v);


// ________________________________________________________________
std::vector<vector2d_t,Eigen::aligned_allocator<vector2d_t>>
    mk_stdvec(const vector2d_t * vb, const vector2d_t * ve);
std::vector<vector3d_t,Eigen::aligned_allocator<vector3d_t>>
    mk_stdvec(const vector3d_t * vb, const vector3d_t * ve);
std::vector<vector4d_t,Eigen::aligned_allocator<vector4d_t>>
    mk_stdvec(const vector4d_t * vb, const vector4d_t * ve);
std::vector<point2d_t,Eigen::aligned_allocator<point2d_t>>
    mk_stdvec(const point2d_t * vb, const point2d_t * ve);

std::vector<point3d_t,Eigen::aligned_allocator<point3d_t>>
    mk_stdvec(const point3d_t * vb, const point3d_t * ve);

std::vector<point4d_t,Eigen::aligned_allocator<point4d_t>>
    mk_stdvec(const point4d_t * vb, const point4d_t * ve);

std::vector<double>
mk_stdvec(const double *vb, const double *ve);
}

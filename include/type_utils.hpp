#ifndef ASTI_TYPE_UTILS
#define ASTI_TYPE_UTILS
#include <type_traits>
#include <Eigen/Core>
template <typename T,int n>
struct alignedType
{
	typedef  T type;
};

template <typename T>
struct alignedType<T,16>
{
	typedef  EIGEN_ALIGN16 T type;
};

#define RAWTYPE(exp)  typename alignedType< typename std::decay<decltype(exp)>::type, __alignof(decltype(exp)) >::type
#endif // ASTI_TYPE_UTILS

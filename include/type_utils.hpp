#ifndef ASTI_TYPE_UTILS
#define ASTI_TYPE_UTILS
#include <type_traits>
#include <Eigen/Core>
#include <vector>
#define STDVEC(TYP)  std::vector<TYP, Eigen::aligned_allocator<TYP>>
#define RAWTYPE(exp) typename std::decay<decltype(exp)>::type
#endif // ASTI_TYPE_UTILS

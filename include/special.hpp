#ifndef ASTI_SPECIAL_HPP
#define ASTI_SPECIAL_HPP
#include <boost/math/special_functions/next.hpp>
namespace util
{
    inline double fprev( double d) { return boost::math::float_prior(d); }
    inline double fnext( double d) { return boost::math::float_next(d); }
}

namespace std
{
#if (_MSC_VER <= 1700)
inline double fmin(double a, double b ) { return a < b? a: b;}
inline double fmax(double a, double b ) { return b < a? a: b;}
#endif
}

#endif // ASTI_SPECIAL_HPP

#ifndef ASTI_SPECIAL_HPP
#define ASTI_SPECIAL_HPP
#include <boost/math/special_functions/next.hpp>
namespace util
{
    inline double fprev( double d) { return boost::math::float_prior(d); }
    inline double fnext( double d) { return boost::math::float_next(d); }
}
#endif // ASTI_SPECIAL_HPP

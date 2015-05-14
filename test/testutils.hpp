#ifndef ASTI_TESTUTILS_HPP
#define ASTI_TESTUTILS_HPP
#include <ostream>
#include "point.hpp"
namespace std {
inline ostream &
operator<<(ostream& os,
           const geom::pt_t<2>& pt)
{
    os << std::setprecision(9) << "(" << pt[0] <<"," << pt[1] << ")";
    return os;
}

inline ostream &
operator<<(ostream& os,
           const geom::vec_t<2>& v)
{
    os << std::setprecision(9) << "(" << v[0] <<"," << v[1] << ")";
    return os;
}
}

#endif // ASTI_TESTUTILS_HPP

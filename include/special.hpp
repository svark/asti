#ifndef ASTI_SPECIAL_HPP
#define ASTI_SPECIAL_HPP
#include <cfloat>

namespace util {
#ifdef _MSC_VER
#if (_MSC_VER <= 1700)
inline double fnext(double a ) { return ::_nextafter(a,a+1) ;}
inline double fprev(double a ) { return ::_nextafter(a,a-1) ;}
#else  //std::nextafter supported in vs 2013
inline double fnext(double a ) { return std::nextafter(a,a+1) ;}
inline double fprev(double a ) { return std::nextafter(a,a-1) ;}
#endif
#else //gcc
inline double fnext(double a ) { return std::nextafter(a,a+1) ;}
inline double fprev(double a ) { return std::nextafter(a,a-1) ;}
#endif
}

namespace std {
#ifdef _MSC_VER // c99 functions not available until vs 12
#if (_MSC_VER <= 1700)
inline double fmin(double a, double b ) { return a < b? a: b;}
inline double fmax(double a, double b ) { return b < a? a: b;}
inline bool isnan(double a)  { return ::_isnan(a) != 0 ;  }
#endif
#endif // _MSC_VER
}

#endif // ASTI_SPECIAL_HPP

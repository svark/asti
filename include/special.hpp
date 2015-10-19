#ifndef ASTI_SPECIAL_HPP
#define ASTI_SPECIAL_HPP
#include <cfloat>
#include "tol.hpp"
#include <limits>

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

inline
double npow(double b, int e)
{
    double result = 1;
    while(e > 0){
        if((e & 1) == 1){
            result *= b;
        }
        b*= b;
        e = e>>1;
    }
    return result;
}

inline
double nroot(double x, long n)
{
    return exp(log(x)/n);
}

}

namespace std {
#ifdef _MSC_VER // c99 functions not available until vs 12
#if (_MSC_VER <= 1700)
inline double fmin(double a, double b ) { return a < b? a: b;}
inline double fmax(double a, double b ) { return b < a? a: b;}
inline bool   isnan(double a)           { return ::_isnan(a) != 0 ;  }
inline double fma(double x, double y, double a) {return x*y+a;}
#endif
#endif // _MSC_VER
}


#endif // ASTI_SPECIAL_HPP

#ifndef ASTI_ARCLENGTH
#define ASTI_ARCLENGTH
#include "integrate1d.hpp"
#include "point.hpp"
#include <functional>
namespace geom {

template <class Curve>
double
arclength(Curve & spl, double t0 = 0.0, double t1 = 1.0)
{
    auto speed = [&spl](double x)->double{
        return len(spl.eval_derivative(1, x));
    };
    return integral(speed, t0, t1);
}

}
#endif // ASTI_ARCLENGTH

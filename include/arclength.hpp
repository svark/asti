#ifndef ASTI_ARCLENGTH
#define ASTI_ARCLENGTH
#include "integrate1d.hpp"
#include "point.hpp"
#include <functional>
namespace geom
{

template <class Curve>
double
arclength(Curve & spl, double t0 = 0.0, double t1 = 1.0)
{
    auto speed = [&spl](double x)->double{
       auto &v = spl.eval_derivative(1, x);
	   return len(v);
    };
    double lastlen = 0;
    double arclen = integrate1d(speed,
                                std::integral_constant < int, 8 >(),
                                t0, t1);
    int n = 8;
    while(tol::neq(lastlen, arclen))
    {
        lastlen = arclen;
        switch(n)
        {
        case 8:
            arclen = integrate1d(speed,
                                 std::integral_constant < int, 16>(),
                                 t0, t1);
            n = 16;
            break;
        case 16:
            arclen = integrate1d(speed,
                                 std::integral_constant < int, 32>(),
                                 t0, t1);
            n = 32;
            break;
        case 32:
            arclen = integrate1d(speed,
                                 std::integral_constant < int, 64>(),
                                 t0, t1);
            n = 64;
            break;
        case 64:
            arclen = integrate1d(speed,
                                 std::integral_constant < int, 128>(),
                                 t0, t1);
            n = 128;
            break;
        case 128:
            arclen = integrate1d(speed,
                                 std::integral_constant < int, 256>(),
                                 t0, t1);
            n = 256;
            break;
        case 256:
            arclen = integrate1d(speed,
                                 std::integral_constant < int, 512>(),
                                 t0, t1);
            n = 512;
            break;
        case 512:
            arclen = integrate1d(speed,
                                 std::integral_constant < int, 1024>(),
                                 t0, t1);
            n = 1024;
            break;
        case 1024:
            break;
        }
    }
    return arclen;
}

}
#endif // ASTI_ARCLENGTH

#ifndef ASTI_INTEGRATE1D_HPP
#define ASTI_INTEGRATE1D_HPP
#include <type_traits>
#include <stdlib.h>
/*
  gl[n_, x_] := Solve[LegendreP[n, x] == 0];
  gldash[n_, a] := D[LegendreP[n, x], x] /. x -> a;
  weights[n_, x_] := 2/((1 - x^2) gldash[n, x]^2);
*/


/*
  Gauss-Legendre n-points quadrature, The method is exact for a
  polynomial of degree <=2n-1 in exact arithmetic
*/

template <int n>
struct gl_points_weights{};

template<>
struct gl_points_weights < 2 > {/* n = 2 */
    static const double points[];
    static const double weights[];
};


template <>
struct gl_points_weights < 4 > {/* n = 4 */
    static const double points[] ;
    static const double weights[] ;
};


template <>
struct gl_points_weights < 8 > {/* n = 8 */
    static const double points[] ;
    static const double weights[] ;
};
template <>
struct gl_points_weights < 16 > {/* n = 16 */
    static const double points[] ;
    static const double weights[] ;
};

template <>
struct gl_points_weights < 32 > {/* n = 32 */
    static const double points[] ;
    static const double weights[] ;
};


template <>
struct gl_points_weights < 64 > {/* n = 64 */
    static const double points[] ;
    static const double weights[] ;
};


template <>
struct gl_points_weights < 128 > {/* n = 128 */ ;
    static const double points[] ;
    static const double weights[] ;
};


template <>
struct gl_points_weights < 256 > {/* n = 256 */
    static const double points[] ;
    static const double weights[] ;
};


template<>
struct  gl_points_weights < 512 > /* n = 512 */ {
    static const double points[] ;
    static const double weights[] ;
};

template <>
struct gl_points_weights < 1024 > {/* n = 1024 */
    static const double points[] ;
    static const double weights[] ;
};


template <class Fn, int n>
double integrate1d(Fn f, std::integral_constant < int, n >,
                   double a,
                   double b)
{
    const double* points = gl_points_weights < n >::points;
    const double* w = gl_points_weights < n >::weights;

    size_t m = n>>1;

    double ga =  0.5 * (b - a);
    double mn =  0.5 * (b + a);

    auto auxf =  [ & f,& ga,& mn](double d) -> double {
        auto gad = ga * d;
        return f(mn + gad) + f(mn - gad);
    };

    double s = 0.0;
    for(size_t i = 0;i < m; ++i)
    {
        s += w[i] * auxf(points[i]);
    }
    return ga*s;
}


template <class Fn>
double integral(Fn speed, double t0, double t1)
{
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

#endif // ASTI_INTEGRATE1D_HPP

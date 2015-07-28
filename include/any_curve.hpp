#ifndef ASTI_ANY_CURVE
#define ASTI_ANY_CURVE
#include <functional>
#include <utility>
#include <type_utils.hpp>
namespace geom
{

template <class Point>
struct  any_curve
{
    template <class Fn>
    any_curve(Fn fn_):fn(fn_)
    {
    }

    Point eval(double u) const
    {
        return fn(u);
    }

    std::function < Point(double u) >  fn;
};

template <class Point>
any_curve < Point> 
make_any_curve(std::function<Point(double u)>&& c)
{
    return any_curve <Point> (std::forward< std::function<Point(double u)> > (c));
}

template <class Curve>
any_curve < RAWTYPE(std::declval<Curve>().eval(0.0)) > 
make_any_curve(Curve c)
{
    typedef RAWTYPE(std::declval<Curve>().eval(0.0))  Point;
    struct Fn
    {
        Fn(Curve c_):c(std::move(c_))
        {
        }
        Point operator()(double u) const
        {
            return c.eval(u);
        }
        Curve c;
    } fn (std::move(c));
    return any_curve <Point> (fn);
}


}
#endif // ASTI_ANY_CURVE

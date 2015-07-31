// -*- mode:c++ -*-
#ifndef ASTI_ANY_CURVE
#define ASTI_ANY_CURVE
#include <functional>
#include <utility>
#include <type_utils.hpp>
namespace geom
{

// 
// abstract curve that uses type erasure to provide an 
// interface consisting of the eval function. Use make_any_curve to turn any
// curve that has eval method into any_curve.
// you could for example provide a quick interface to the twisted cubic (t, t^2,
// t^3) by doing the following
// std::function < point3d_t (double) > fn = [](double t){return
// make_pt(t,t * t, t * t * t);}
// auto ac = make_any_curve(fn)
//
template <class Point>
struct  any_curve
{
    template <class Fn>
    any_curve(Fn fn_):fn(fn_)
    {
    }
	any_curve(const any_curve& o)
		:  fn(o.fn){}

	any_curve(any_curve&& o):
		fn(std::forward< std::function<Point (double) > >( o.fn) )
	{
	}
    Point eval(double u) const
    {
        return fn(u);
    }

    std::function < Point(double ) >  fn;
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

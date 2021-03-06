#ifndef ASTI_UTIL_HPP
#define ASTI_UTIL_HPP
#include <vector>
#include <iterator>
#include <functional>
#include <utility>

namespace util {

enum infix1 { inside , outside, boundary_of};

template <class Fn>
struct infix_op
{
    infix_op(Fn f_):f(f_)
    {
    }

    Fn f;
};

template <class Scalar>
struct curried_arg_in
{
    curried_arg_in(Scalar x_):x(x_){}
    bool operator()( const std::pair<Scalar,Scalar>& y) const
    {
        return x > y.first  && x < y.second ;
    }
    Scalar x;
};
template <class Scalar>
struct curried_arg_out
{
    curried_arg_out(Scalar x_):x(x_){}
    bool operator()( const std::pair<double,double>& y) const
    {
        return x < y.first  || x > y.second ;
    }
    Scalar x;
};
template <class Scalar>
struct curried_arg_bndry
{
    curried_arg_bndry(Scalar x_):x(x_){}
    bool operator()( const std::pair<Scalar,Scalar>& y) const
    {
        return x == y.first  || x == y.second ;
    }
    Scalar x;
};

template <typename Scalar>
infix_op <curried_arg_in<Scalar> >
operator < (Scalar x, std::integral_constant<infix1,inside> i)
{
    return infix_op<curried_arg_in<Scalar>> (curried_arg_in<Scalar>(x));
}
template <typename Scalar>
infix_op <curried_arg_out<Scalar> >
operator < (Scalar x, std::integral_constant<infix1,outside> i)
{
    return infix_op<curried_arg_out<Scalar> > (curried_arg_out<Scalar> (x));
}

template <typename Scalar>
infix_op <curried_arg_bndry<Scalar> >
operator < (Scalar x, std::integral_constant<infix1,boundary_of> i)
{
    return infix_op<curried_arg_bndry<Scalar>> (curried_arg_bndry<Scalar>(x));
}

template <class Fn,class Scalar>
bool operator > ( infix_op<Fn> op,
                  const std::pair<Scalar,Scalar>&  y)
{
    return op.f(y);
}



template<class Fn,class A, class AllocT>
std::vector<decltype(Fn()(A()))> fmap(Fn f, const std::vector<A, AllocT >& as) {
    typedef decltype(Fn()(A())) B;
    std::vector<B> result;
    std::transform(as.begin(), as.end(),
                   std::back_inserter(result), f);
    return result;
}

struct tag_switcher
{
    template <class Fn1, class Fn2>
    static auto eval(Fn1 f1_, Fn2 f2_, std::false_type) { return f2_();}

    template <class Fn1, class Fn2>
    static auto eval(Fn1 f1_, Fn2 f2_, std::true_type)  { return f1_();}

};

using std::next;
}


static const std::integral_constant<util::infix1,util::boundary_of> _on_ =  std::integral_constant<util::infix1,util::boundary_of>();
static const std::integral_constant<util::infix1,util::inside>      _in_ = std::integral_constant<util::infix1,util::inside>();
static const std::integral_constant<util::infix1,util::outside>     _out_= std::integral_constant<util::infix1,util::outside>();

#endif//UTIL_HPP

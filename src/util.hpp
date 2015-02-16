//#include <boost/util.hpp>
#include <list>
#include <vector>
#include <iterator>
#include <functional>

namespace util {

#ifdef VARIADIC_ARGS_SUPPORTED
template <class FnType,class Args...>
void switch_bool(b,f,args...) {
    if(b)
        f(args, boost::mpl::integral_c<true>());
    else
        f(args, boost::mpl::integral_c<false>());
}
#else
#define switch_bool(b,f,...) {                                      \
        if(b)                                                       \
            f(__VA_ARGS__, boost::mpl::integral_c<bool,true>());    \
        else                                                        \
            f(__VA_ARGS__, boost::mpl::integral_c<bool,false>());   \
    }

#define switch_case_3_inner(EnumT,e,lo,f,...)                   \
    case lo:                                                    \
    f(__VA_ARGS__, boost::mpl::integral_c<EnumT,EnumT(lo)>());  \
    break;                                                      \
case (lo+1):                                                    \
f(__VA_ARGS__, boost::mpl::integral_c<EnumT,EnumT(lo+1)>());    \
break;                                                          \
case lo+2:                                                      \
f(__VA_ARGS__, boost::mpl::integral_c<EnumT,EnumT(lo+2)>());    \
break;

#define switch_case_4_inner(EnumT,e,lo,f,...)                       \
    switch_case_3_inner(EnumT,e,lo,f,__VA_ARGS__)                   \
case lo+3:                                                          \
    f(__VA_ARGS__, boost::mpl::integral_c<EnumT,EnumT(lo+3)>());    \
    break;                                                          \

#define switch_case_5_inner(EnumT,e,lo,f,...)                       \
    switch_case_4_inner(EnumT,e,lo,f,__VA_ARGS__)                   \
case lo+4:                                                          \
    f(__VA_ARGS__, boost::mpl::integral_c<EnumT,EnumT(lo+3)>());    \
    break;                                                          \

#define switch_case_3(EnumT,e,lo,err,f,...)  /* */      \
    switch(e)                                           \
    {                                                   \
        switch_case_3_inner(EnumT,e,lo,f,__VA_ARGS__)   \
    default:                                            \
            err();                                      \
    }
#define switch_case_4(EnumT,e,lo,err,f,...)             \
    switch(e)                                           \
    {                                                   \
        switch_case_4_inner(EnumT,e,lo,f,__VA_ARGS__)   \
    default:                                            \
            err();                                      \
    }

#define switch_case_5(EnumT,e,lo,err,f,...)             \
    switch(e)                                           \
    {                                                   \
        switch_case_5_inner(EnumT,e,lo,f,__VA_ARGS__)   \
    default:                                            \
            err();                                      \
    }
#endif

// utility class to allow range based for
template <class T>
struct rvector
{
    typedef std::reverse_iterator<typename std::vector<T>::const_iterator > CIT;
    typedef std::reverse_iterator<typename std::vector<T>::iterator > IT;
    rvector(std::vector<T>& v_):v(v_){}
    T operator[](int i) { return v.rbegin()[i];}
    CIT cbegin() const { return v.rcbegin(); }
    CIT cend() const { return v.rcbegin(); }
    IT begin() { return v.rbegin(); }
    IT end()   { return v.rend(); }
private:
    std::vector<T> & v;
};

template <class IterT>
struct range
{
    range(IterT b_, IterT e_):b(b_), e(e_){}
    typename std::iterator_traits < IterT >::value_type operator[](int i) { return b[i];}
    IterT begin() { return b; }
    IterT end() { return e; }
private:
    IterT b, e;
};

template<class IterT>
range < IterT > mk_range(IterT b, IterT e) { return range<IterT>(b, e);}


template <class IterT>
struct rrange
{
    rrange(IterT b_, IterT e_):b(e_), e(b_){}
    typename std::iterator_traits < IterT >::value_type  operator[](int i) { return b[i];}
    std::reverse_iterator<IterT> begin() { return b; }
    std::reverse_iterator<IterT> end() { return e; }
private:
    std::reverse_iterator<IterT> b, e;
};

template<class IterT>
rrange < IterT > mk_rrange(IterT b, IterT e) { return rrange<IterT>(b, e);}

template <class knots_t>
struct skip_ith_iter
{
    skip_ith_iter(int i_, knots_t& knots_)
        :knots(knots_), i(i_){

    }

    skip_ith_iter &operator++() {
        ++iter;
        if( std::distance(knots.begin(), iter) == i )
            ++iter;
        return *this;
    }

    double operator *() { return *iter;}
private:
    const knots_t& knots;
    int i;
    typename knots_t::const_iterator iter;
};

template <class ContT>
struct back_n_insert_iterator : std::back_insert_iterator<ContT>
{
    back_n_insert_iterator(int n, ContT&v):
        std::back_insert_iterator<ContT>(v),_n(n)
    {
    }
    typedef typename ContT::value_type vt_t;
    back_n_insert_iterator& operator=(const vt_t& v)
    {
        for(int i =0; i < _n;++i)
            container->push_back(v);
        return *this;
    }
    back_n_insert_iterator& operator*()
    {   // pretend to return designated value
        return (*this);
    }

    back_n_insert_iterator& operator++()
    {   // pretend to preincrement
        return (*this);
    }

    back_n_insert_iterator operator++(int)
    {   // pretend to postincrement
        return (*this);
    }
    int _n;
};

template <class ContT>
back_n_insert_iterator<ContT> back_n_inserter(int n, ContT&v) { 
    return (back_n_insert_iterator<ContT>(n,v)); 
}

template<class A, class Fn>
std::vector<decltype(Fn()(A()))> fmap(Fn f, std::vector<A> as) {
    typedef decltype(Fn()(A())) B;
    std::vector<B> result;
    std::transform(as.begin(), as.end(),
                   std::back_inserter(result), f);
    return result;
}
template<class A, class Fn>
std::list<decltype(Fn()(A()))> fmap(Fn f, std::list<A> as ) {
    std::list<decltype(Fn()(A()))> result;
    std::transform(as.begin(), as.end(),
                   std::back_inserter(result), f);
    return result;
}

using std::next;
}
#include <utility>
template <class ContT>
struct std::iterator_traits<util::back_n_insert_iterator<ContT>> :
    std::iterator_traits<std::back_insert_iterator<ContT> >
{
};
#ifndef _MSC_VER
namespace stdext{
template<class ArrayT>
ArrayT * make_checked_array_iterator(ArrayT * arr, int size) { return arr;}
}
#endif

namespace std {
template<class _Container>
struct _Is_checked_helper<util::back_n_insert_iterator<_Container> >
    : public true_type
{   // mark back_insert_n_iterator as checked
};
}

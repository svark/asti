//#include <boost/util.hpp>
#include <list>
#include <vector>
#include <iterator>
#include <functional>
namespace util {

#ifdef VARIADIC_ARGS_SUPPORTED
  template <class FnType,class Args...>
  void switch_on(b,f,args...) {
    if(b)
      f(args, boost::mpl::integral_c<true>());
    else
      f(args, boost::mpl::integral_c<false>());
  }
#else
#define switch_bool(b,f,...) {                                    \
        if(b)                                                   \
            f(__VA_ARGS__, boost::mpl::integral_c<bool,true>());     \
        else                                                    \
            f(__VA_ARGS__, boost::mpl::integral_c<bool,false>());    \
  }

#define switch_case_3_inner(EnumT,e,lo,f,...)\
    case lo:                                                  \
        f(__VA_ARGS__, boost::mpl::integral_c<EnumT,EnumT(lo)>());  \
        break;                                                \
    case (lo+1):                                              \
        f(__VA_ARGS__, boost::mpl::integral_c<EnumT,EnumT(lo+1)>());\
        break;                                                \
    case lo+2:                                                \
        f(__VA_ARGS__, boost::mpl::integral_c<EnumT,EnumT(lo+2)>());\
        break;

#define switch_case_4_inner(EnumT,e,lo,f,...)\
 switch_case_3_inner(EnumT,e,lo,f,__VA_ARGS__)                 \
    case lo+3:                                                 \
        f(__VA_ARGS__, boost::mpl::integral_c<EnumT,EnumT(lo+3)>()); \
        break;                                                 \

#define switch_case_5_inner(EnumT,e,lo,f,...)\
 switch_case_4_inner(EnumT,e,lo,f,__VA_ARGS__)                 \
    case lo+4:                                                 \
        f(__VA_ARGS__, boost::mpl::integral_c<EnumT,EnumT(lo+3)>()); \
        break;                                                 \


#define switch_case_3(EnumT,e,lo,err,f,...)  /* */            \
    switch(e)                                                 \
    {                                                         \
        switch_case_3_inner(EnumT,e,lo,f,__VA_ARGS__)         \
    default:                                                  \
        err();                                                \
    }
#define switch_case_4(EnumT,e,lo,err,f,...)                   \
    switch(e)                                                 \
    {                                                         \
        switch_case_4_inner(EnumT,e,lo,f,__VA_ARGS__)         \
    default:                                                  \
        err();                                                \
    }

#define switch_case_5(EnumT,e,lo,err,f,...)                   \
    switch(e)                                                 \
    {                                                         \
        switch_case_5_inner(EnumT,e,lo,f,__VA_ARGS__)         \
    default:                                                  \
        err();                                                \
    }
#endif

    template <class T>
    struct rvector
    {
        rvector(std::vector<T>& v_):v(v_){}
        T operator[](int i) { return v.rbegin()[i];}
        std::vector<T> & v;
    };
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

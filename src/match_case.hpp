#ifndef MATCH_CASE_HPP
#define MATCH_CASE_HPP
namespace util
{
    #ifdef VARIADIC_ARGS_SUPPORTED
template <class FnType,class Args...>
void switch_bool(b,f,args...) {
    if(b)
        f(args, std::integral_constant<true>());
    else
        f(args, std::integral_constant<false>());
}
#else
#define switch_bool(b,f,...) {                                      \
        if(b)                                                       \
            f(__VA_ARGS__, std::integral_constant<bool,true>());    \
        else                                                        \
            f(__VA_ARGS__, std::integral_constant<bool,false>());   \
    }

#define switch_bool_ret(b,f,...) {                                      \
        if(b)                                                       \
            return f(__VA_ARGS__, std::integral_constant<bool,true>());    \
        else                                                        \
            return f(__VA_ARGS__, std::integral_constant<bool,false>());   \
    }

#define switch_case_3_inner(EnumT,e,lo,f,...)                   \
    case lo:                                                    \
    f(__VA_ARGS__, std::integral_constant<EnumT,EnumT(lo)>());  \
    break;                                                      \
case (lo+1):                                                    \
f(__VA_ARGS__, std::integral_constant<EnumT,EnumT(lo+1)>());    \
break;                                                          \
case lo+2:                                                      \
f(__VA_ARGS__, std::integral_constant<EnumT,EnumT(lo+2)>());    \
break;

#define switch_case_4_inner(EnumT,e,lo,f,...)                       \
    switch_case_3_inner(EnumT,e,lo,f,__VA_ARGS__)                   \
case lo+3:                                                          \
    f(__VA_ARGS__, std::integral_constant<EnumT,EnumT(lo+3)>());    \
    break;                                                          \

#define switch_case_5_inner(EnumT,e,lo,f,...)                       \
    switch_case_4_inner(EnumT,e,lo,f,__VA_ARGS__)                   \
case lo+4:                                                          \
    f(__VA_ARGS__, std::integral_constant<EnumT,EnumT(lo+3)>());    \
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

}
#endif

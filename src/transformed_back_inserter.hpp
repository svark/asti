#ifndef TRANSFORMED_BACK_INSERTER_HPP
#define TRANSFORMED_BACK_INSERTER_HPP
#include <iterator>
namespace util {
template <class ContT, class Ret, class Arg>
struct transformed_back_insert_iterator : std::back_insert_iterator<ContT>
{
    transformed_back_insert_iterator(ContT&v, Ret ( * f)(const Arg & )):
        std::back_insert_iterator<ContT>(v),_f(f)
    {
    }
    transformed_back_insert_iterator& operator=(const Arg& v)
    {
        container->push_back(_f(v));
        return *this;
    }
    transformed_back_insert_iterator& operator*()
    {   // pretend to return designated value
        return (*this);
    }

    transformed_back_insert_iterator& operator++()
    {   // pretend to preincrement
        return (*this);
    }

    transformed_back_insert_iterator operator++(int)
    {   // pretend to postincrement
        return (*this);
    }
    Ret ( * _f)(const Arg & );
};

template <class ContT, class Ret, class Arg>
transformed_back_insert_iterator<ContT,Ret,Arg>
transformed_back_inserter(ContT&v,
                          Ret (*f)(const Arg &))
{
    return (transformed_back_insert_iterator<ContT, Ret, Arg>(v, f));
}

}

namespace std {
template <class ContT, class Ret, class Arg>
struct iterator_traits< util::transformed_back_insert_iterator<ContT, Ret, Arg> >:
        std::iterator_traits<std::back_insert_iterator<ContT> >
{
};

template<class _Container, class Ret, class Arg>
struct _Is_checked_helper<util::transformed_back_insert_iterator<_Container,
                                                                 Ret, Arg> >
    : public true_type
{   // mark back_insert_n_iterator as checked
};

}

#endif //TRANSFORMED_BACK_INSERTER_HPP

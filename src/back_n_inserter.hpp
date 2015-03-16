#ifndef BACK_N_INSERTER_HPP
#define BACK_N_INSERTER_HPP
namespace util {
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


template <class ContT>
struct std::iterator_traits<util::back_n_insert_iterator<ContT>> :
    std::iterator_traits<std::back_insert_iterator<ContT> >
{
};
}



namespace std {
template<class _Container>
struct _Is_checked_helper<util::back_n_insert_iterator<_Container> >
    : public true_type
{   // mark back_insert_n_iterator as checked
};
}

#endif //BACK_N_INSERTER_HPP

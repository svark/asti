#ifndef ASTI_RANGE_HPP
#define ASTI_RANGE_HPP
#include <iterator>
namespace util
{


// utility class to allow range based for
template <class VecT>
struct rvector
{
    typedef std::reverse_iterator<typename VecT::const_iterator > CIT;
    typedef std::reverse_iterator<typename VecT::iterator > IT;
    typedef typename VecT::value_type T;
    rvector(VecT& v_):v(v_){}
    T operator[](int i) { return v.rbegin()[i];}
    CIT cbegin() const { return v.rcbegin(); }
    CIT cend() const { return v.rcbegin(); }
    IT begin() { return v.rbegin(); }
    IT end()   { return v.rend(); }
private:
    VecT & v;
};


template <class IterT>
struct range
{
    range(IterT b_, IterT e_):b(b_), e(e_){}
    typename std::iterator_traits < IterT >::value_type
    operator[](int i) {return b[i];}
    IterT begin() const { return b; }
    IterT end() const { return e; }

    IterT cbegin() const { return b; }
    IterT cend() const { return e; }

    std::reverse_iterator<IterT> rbegin() const {
        return std::reverse_iterator<IterT>(e);
    }
    std::reverse_iterator<IterT> rend() const {
        return std::reverse_iterator<IterT>(b);
    }

    typename std::iterator_traits < IterT >::value_type
    front() const { return b[0] ; }

    typename std::iterator_traits < IterT >::value_type
    back() const { return e[-1] ; }

    size_t size() const { return std::distance(b, e);}
private:
    IterT b, e;
};

template<class IterT>
range < IterT > mk_range(IterT b, IterT e) { return range<IterT>(b, e);}


template <class IterT>
struct rrange
{
    rrange(IterT b_, IterT e_):b(e_), e(b_){}
    typename std::iterator_traits < IterT >::value_type
    operator[] (int i) {return b[i];}
    std::reverse_iterator<IterT> begin() { return b; }
    std::reverse_iterator<IterT> end() { return e; }
private:
    std::reverse_iterator<IterT> b, e;
};

template<class IterT>
rrange < IterT > mk_rrange(IterT b, IterT e) { return rrange<IterT>(b, e);}

}

#endif//RANGE_HPP

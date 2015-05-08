#ifndef ASTI_CONSTANT_ITER_HPP
#define ASTI_CONSTANT_ITER_HPP
#include <iterator>
namespace util {
template <typename T>
struct constant_iterator:
        std::iterator<std::input_iterator_tag, T>
{
	constant_iterator(T v_,std::ptrdiff_t idx_=0):
        v(v_), idx(idx_)
    {
    }
	constant_iterator(const constant_iterator&other)
		:v(other.v),idx(other.idx)
	{
	}

    T operator*() const
    {
        return v;
    }

	constant_iterator& operator+=(std::ptrdiff_t sz)
    {
		idx+=sz;
        return *this;
    }

    T operator[](size_t i)  const { return v;}

    constant_iterator& operator++()
    {   // pretend to preincrement
		++idx;	
        return (*this);
    }
	
    constant_iterator operator++(int)
    {   // pretend to postincrement
		constant_iterator tmp(*this);
		++idx;
		return tmp;
    }
	constant_iterator& operator--()
    {   // pretend to preincrement
		--idx;	
        return (*this);
    }
	
    constant_iterator operator--(int)
    {   // pretend to postincrement
		constant_iterator tmp(*this);
		--idx;
		return tmp;
    }
    T v;
	std::ptrdiff_t idx;
};

template <class T>
constant_iterator<T> operator+(constant_iterator<T>& it, std::ptrdiff_t sz) { 
	constant_iterator<T> it_(it);
	it_+=sz;
	return it_; 
}

template <class T>
bool
operator==(const constant_iterator<T>&it1, const constant_iterator<T>&it2)
{
	return std::make_pair(it1.v,it1.idx)==std::make_pair(it2.v,it2.idx);
}
template <class T>
bool
operator!=(const constant_iterator<T>&it1, const constant_iterator<T>&it2)
{
	return std::make_pair(it1.v,it1.idx)!=std::make_pair(it2.v,it2.idx);
}

template <class T>
constant_iterator<T>  make_constant_iterator( T v, std::ptrdiff_t idx=0) { 
	return constant_iterator<T>(v, idx) ; 
}
}
#endif //CONSTANT_ITER_HPP

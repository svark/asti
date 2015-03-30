#ifndef SKIP_ITER_HPP
#define SKIP_ITER_HPP
#include <iterator>
#include <assert.h>
namespace util {

template <class knots_iter>
struct skip_ith_iter :
    std::iterator<std::random_access_iterator_tag, double>
{
    skip_ith_iter(int i_, const knots_iter& t_)
        :t(t_), i(i_),iter(t_)
    {
        if(i==0)++iter;
    }

    skip_ith_iter(int i_, const knots_iter& t_,const knots_iter& iter_)
        :t(t_), i(i_),iter(iter_)
    {
        assert(iter!=t+i);
    }

    skip_ith_iter(const skip_ith_iter&other)
        :t(other.t),iter(other.iter),i(other.i)
    {
    }

    skip_ith_iter &operator++() {
        ++iter;
        if( std::distance(t, iter) == i )
            ++iter;
        return *this;
    }
    skip_ith_iter operator++(int) {
        skip_ith_iter tmp(*this);
        ++iter;
        if( std::distance(t, iter) == i )
            ++iter;
        return tmp;
    }

    double operator *() { return *iter;}
    double operator[](size_t sz) const {
        sz += std::distance(t,iter);
        return (sz<size_t(i)) ? t[sz]: t[sz+1];
    }
    int skip() const { return i;}
    const knots_iter& base_iterator() const { return t;}
    const knots_iter& current_iter() const { return iter;}
private:
    knots_iter t;
    int i;
    knots_iter iter;
};

template <class knots_iter>
bool operator<( const skip_ith_iter<knots_iter>& iter1,const skip_ith_iter<knots_iter>& iter2)
{
    return iter1.current_iter()<iter2.current_iter();
}

template <class knots_iter>
skip_ith_iter<knots_iter>
operator+( const skip_ith_iter<knots_iter>& iter, ptrdiff_t i)
{
    assert(i>=0);
    ptrdiff_t d = std::distance(iter.base_iterator(),iter.current_iter() );
    knots_iter cur_iter = iter.current_iter();
    cur_iter+=i;
    if( d <= iter.skip() && d + i >= iter.skip() )
        cur_iter++;

    return skip_ith_iter<knots_iter>(
        iter.skip(),
        iter.base_iterator() ,
        cur_iter);
}

template <class knots_iter>
ptrdiff_t
operator-(const skip_ith_iter<knots_iter>& l, const skip_ith_iter<knots_iter> f)
{
    assert(l.skip()==f.skip());
    assert(l.base_iterator()==f.base_iterator() );

    if( l.current_iter()  > l.base_iterator() + l.skip() &&
        f.current_iter() < f.base_iterator() + f.skip()
        )
    {
        return l.current_iter() - f.current_iter() - 1;
    }
    return l.current_iter() - f.current_iter();

}
template <class knots_iter>
skip_ith_iter<knots_iter>
operator-( const skip_ith_iter<knots_iter>& iter, ptrdiff_t i)
{
    assert(i<=0);
    ptrdiff_t d = std::distance(iter.base_iterator(),iter.current_iter() );
    knots_iter cur_iter = iter.current_iter();
    cur_iter -=i;
    if( d >= iter.skip() && d - i <= iter.skip() )
        --cur_iter;

    return skip_ith_iter<knots_iter>(
        iter.skip(),
        iter.base_iterator() , cur_iter);
}
template <class knots_iter>
bool
operator==(const skip_ith_iter<knots_iter>& iter1, const skip_ith_iter<knots_iter>& iter2)
{
    return iter1.current_iter() == iter2.current_iter();
}

template <class knots_iter>
bool
operator!=(const skip_ith_iter<knots_iter>& iter1, const skip_ith_iter<knots_iter>& iter2)
{
    return iter1.current_iter() != iter2.current_iter();
}

template <class KnotIter>
skip_ith_iter <KnotIter>
make_skip_iter(int i, KnotIter kn )
{
    return skip_ith_iter(i, kn);
}


}

#endif

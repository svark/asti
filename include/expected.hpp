#ifndef ASTI_EXPECTED
#define ASTI_EXPECTED
#include "geom_exception.hpp"
#include <assert.h>
#include "type_utils.hpp"
namespace util
{
using geom::geom_error_code_t;
struct unexpected
{
    geom_error_code_t e;
};

inline unexpected
make_unexpected(geom_error_code_t e)
{
    unexpected unex = {e};
    return unex;
}


template <class U>
class expected
{
    typename std::aligned_union<sizeof(U), U, geom_error_code_t>::type  storage;
    bool valid_;
public:
    typedef U value_type;
    expected(const U &t)
    {
        new (&storage) U(t);
        valid_ = true;
    }

    expected(U &&t)
    {
        new (&storage) U(std::forward<U>(t));
        valid_ = true;
    }

    expected(const expected& o)
    {
        if(!!o)  {
            new (&storage) U(*o);
            valid_ = true;
        }
        else {
            new (&storage) geom_error_code_t(o.code());
            valid_ = false;
        }

    }
    expected(expected&& o)
    {
        valid_ = o.valid();
        if(!!o)
            new (&storage) U(std::move(*o));
        else
            new (&storage) geom_error_code_t(o.code());
    }
    expected & operator=(expected& other)
    {
        if(valid())
        {
            (**this).~U();
        }
        valid_ = o.valid();
        if(!!o)
            new (&storage) U(*o);
        else
            new (&storage) geom_error_code_t(o.code());
        return *this;
    }
    expected(const unexpected& unex)
    {
        new (&storage) geom_error_code_t(unex.e);
        valid_ = false;
    }

    U& operator*()  {
        assert(valid_);
        return  *reinterpret_cast<U*>(&storage);
    }

    const U& operator*() const {
        assert(valid_); return *reinterpret_cast<const U*>(&storage);
    }

    operator bool() const { return valid_; }
    void valid() const { return valid_; }

    geom_error_code_t code() const {
        assert(!valid_);
        return *reinterpret_cast<const geom_error_code_t* >(&storage);
    }
    ~expected()
    {
        if(valid_)
            reinterpret_cast<U*>(&storage)->~U();
    }
};

template <>
class expected<void>
{
    geom_error_code_t e_;
public:
    typedef void value_type;
    expected():e_(geom::no_error)
    {

    }
    expected(const unexpected& unex) : e_(unex.e)
    {
    }

    void operator*() const {    }
    operator bool() const { return valid(); }
    bool valid() const { return e_==0; }
    geom_error_code_t code() const { assert(!valid()); return e_;}
};

#define ASTI_EXPECT_RET(ex,expr)                \
    if(!!ex) {                                  \
        return expr;                            \
    }                                           \
    return util::make_unexpected(ex.code());

#define ASTI_EXPECT_DO(ex,expr)                 \
    if(!!ex) {                                  \
        expr;                                   \
    }                                           \
    return util::make_unexpected(ex.code());

#define ASTI_CHECK(ex)                                  \
    if(!ex) return util::make_unexpected(ex.code());

static const expected<void> exvoid;

template <class ExpectedT,class Fn>
auto bind_expected(ExpectedT ex, Fn f) -> expected<RAWTYPE(f(*ex))>
{
    if(!!ex) {
        return f(*ex);
    }
    return make_unexpected(ex.code());
};


template <class Fn>
auto bind_expected(expected<void> ex, Fn f) -> expected<RAWTYPE(f())>
{
    if(!!ex) {
        return f();
    }
    return make_unexpected(ex.code());
};
}
#endif // ASTI_EXPECTED

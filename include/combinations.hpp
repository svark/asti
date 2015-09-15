#ifndef ASTI_COMBINATIONS_HPP
#define ASTI_COMBINATIONS_HPP

namespace util {


// return array of C^n_k for k in [0,n]
inline
std::vector<size_t> ncks(size_t n)
{
    std::vector<size_t> ncks_;
    ncks_.reserve(n + 1);
    size_t val = 1;
    ncks_.push_back(1);
    for(size_t k = 1; k <= n; ++k)
    {
        val *= (n - k + 1);
        val /= k;
        ncks_.push_back(val);
    }
    return ncks_;
}

// return C^{n+k}_k values for k in [0,j]
inline
std::vector<size_t> nkcks(size_t n, size_t j)
{
    std::vector<size_t> nkcks_;
    nkcks_.reserve(j + 1);
    size_t val = 1;
    nkcks_.push_back(1);
    for(size_t k = 1; k <= j; ++k)
    {
        val *= (n + k);
        val /= k;
        nkcks_.push_back(val);
    }
    return nkcks_;
}
}
#endif // ASTI_COMBINATIONS_HPP

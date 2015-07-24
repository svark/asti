#ifndef ASTI_RMAT_EXPLICIT
#define ASTI_RMAT_EXPLICIT
#include "tol.hpp"
#include "type_utils.hpp"
namespace geom {
// -- knot insertion helper classes
// computes rmat_explicit entries for the oslo algo.
// 1< k < d that identifies the R-matrix (&its size k * k+1)
// x new value of knot to insert
// mu - index to knot vector so that t[u] <= x < t[u+1]
// i , j <-- matrix entry we are interested in
// for the oslo algorithm & and also for finding derivatives

inline static double sdiv(double n, double d) {
    if (tol::param_eq(d, 0) && tol::param_eq(n, 0))
        return 0.0;
    else
        return n / d;
}
inline static double sdiv1(double n, double d) {
    if (tol::param_eq(d, 0) && tol::param_eq(n, 0))
        return 1.0;
    else
        return n / d;
}

template <class KnotIter>
struct rmat_explicit
{
    rmat_explicit(
        KnotIter knots_, // knot iterator after which
                         // to insert
        int      d,      // 1  <= k <= degree
        double   u       // knot to insert
        ):t(knots_), k(d), tau(u)
    {
    }

    static double sdiv(double n, double d) {
        if(tol::param_eq(d, 0))
            return 0;
        return n / d;
    }

    double get_diag(size_t i) const
    {
        double denom = (t[i + 1] - t[i + 1 - k]);
        return sdiv(t[i + 1] - tau, denom);
    }

    double get_ndiag(size_t i) const
    {
        double denom = (t[i + 1] - t[i + 1 - k]);
        return sdiv(tau - t[i + 1 - k], denom);
    }

    inline int size() const
    {
        return k + 1;
    }

    KnotIter t;                 // knots;
    int      k;                 // size
    double   tau;               // knot to insert

private:
    rmat_explicit & operator=(const rmat_explicit & other);
};


// derivatives.
template <class KnotIter>
class der_rmat_explicit:
        public rmat_explicit<KnotIter>
{
    using rmat_explicit<KnotIter>::t;
    using rmat_explicit<KnotIter>::k;
    using rmat_explicit<KnotIter>::tau;
public:
    der_rmat_explicit(
        KnotIter knots_iter,    // knot
        // iterator
        // after which
        // to insert
        int d ,                // degree
        double u             // knot to insert
        ):rmat_explicit<KnotIter>(knots_iter, d, u)
    {
    }

    double get_diag(size_t i) const
    {
        double denom = (t[i+1] - t[i+1 - k]);
        return -1/ denom;
    }
    double get_ndiag(size_t i) const
    {
        double denom = (t[i+1] - t[i+1 - k]);
        return 1/ denom;
    }
};

// mult_rmat provides temp storage for multiplying knot insertion
// matrices since multiplication happens left to right its row size is
// always 1
struct mult_rmat
{
    mult_rmat():basis(1, 1)      // default is the identity matrix of
                                 // dim 1
    {
    }

public:
    template < class KnotInsertionMatrix >
    mult_rmat(const mult_rmat & m,
              const KnotInsertionMatrix & kim):basis(kim.size())
    {
        size_t num_cols = kim.size();
        assert(m.basis.size() + 1 == num_cols); // mult order
        basis[0] = (m.basis[0] * kim.get_diag(0));
        for (size_t k = 1; k < num_cols - 1; ++k) {
            basis[k] = (m.basis[k - 1] * kim.get_ndiag(k - 1))
                +(m.basis[k]  * kim.get_diag(k));
        }
        basis[num_cols - 1] =
            m.basis[num_cols - 2] * kim.get_ndiag(num_cols - 2);
    }

    mult_rmat(mult_rmat && other):basis(std::forward <std::vector<double>>
                                        (other.basis))
    {
    }
    void swap(mult_rmat & m)
    {
        basis.swap(m.basis);
    }

    template < class KnotInsertionMatrix >
    mult_rmat & operator*=(const KnotInsertionMatrix & kim)
    {
        mult_rmat m(*this, kim);
        m.swap(*this);
        return *this;
    }

    double get(size_t k) const
    {
        return basis[k];
    }

    size_t size() const
    {
        return basis.size();
    }

    const std::vector<double>&
    getb() const {
        return basis;
    }

    std::vector<double>&&
    moveb() {
        return std::move(basis);
    }

    template<class PointIter>
    auto prod(const PointIter cpts) ->
        RAWTYPE(cpts[0])
    {
        typedef  RAWTYPE(cpts[0]) point_type;
        point_type pt(0.0);
        for (int k = 0; k < size(); ++k)
        {
            pt += get(k) * cpts[k];
        }
        return pt;
    }
private:
    std::vector < double > basis;
};

template <class KnotIter>
rmat_explicit < KnotIter >
make_rmat_explicit(KnotIter t, int d, double u)
{
    return rmat_explicit < KnotIter > (t, d, u);
}

template <class KnotIter>
rmat_explicit < KnotIter >
make_der_rmat_explicit(KnotIter t, int d, double u)
{
    return der_rmat_explicit < KnotIter > (t, d, u);
}


} // namespace geom
#endif // ASTI_RMAT_EXPLICIT

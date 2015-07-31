//-*- mode:c++ -*-
#ifndef ASTI_MONOMIAL_FORM
#define ASTI_MONOMIAL_FORM
#include "type_utils.hpp"
#include "point.hpp"
namespace geom {

template <class Point>
struct monomial_form
{
    typedef RAWTYPE(mk_stdvec(Point(0.0))) cpts_t;
    monomial_form(cpts_t a_):a(std::move(a_)){}

    Point eval(double u) {
        assert(a.size());
        Point v(0.0);
        double ui = 1;
        for(int i = 0;i < a.size(); ++i, ui *= u)
        {
            v += ui * a[i];
        }
        return v;
    }
    Point operator[](int i) const { return a[i];}
    size_t size() const { return a.size(); }
    int    degree() const { return int(a.size() - 1);}
    cpts_t const & coeffs() const { return a;}

private:
    cpts_t a;
};
}
#endif // ASTI_MONOMIAL_FORM

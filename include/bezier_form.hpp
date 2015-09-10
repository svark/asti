#ifndef ASTI_BEZIER_FORM
#define ASTI_BEZIER_FORM
#include <vector>
#include "point_fwd.hpp"
#include "bspline.hpp"
#include "bspline_queries.hpp"

namespace geom {

// though the name says "bezier form" it is actually a bspline,  the key
// difference being that outside the iterval of definition the values
// returned by eval function are zero.
template <class Point>
struct bezier_form : public bspline < Point >
{
    typedef bspline<Point> base_t;
    bezier_form(RAWTYPE(mk_stdvec(Point(0.0))) c,
                double a = 0,  double b =  1.0):
        base_t(std::move(c),
               std::vector<double>(2 * c.size()),
               c.size() - 1)
    {
        int deg =  base_t::degree();
        std::fill_n(base_t::t.begin() + (      0), deg + 1, a);
        std::fill_n(base_t::t.begin() + (deg + 1), deg + 1, b);
    }

    bezier_form(base_t && b):base_t(std::forward < base_t > (b))
    {
        check_invariants();
    }
    void check_invariants()
    {
        base_t::check_invariants();
        assert(qry::is_bezier(*this));
    }
};

template <class CptsT>
auto make_bezier_form(CptsT cpts,
                      double a,  double b) -> bezier_form<RAWTYPE(cpts[0])>
{
    return bezier_form<RAWTYPE(cpts[0])>(std::move(cpts), a, b);
}

}

#endif // ASTI_BEZIER_FORM

#ifndef ASTI_BEZIER_FORM
#define ASTI_BEZIER_FORM
#include <vector>
#include "point_fwd.hpp"
#include "bspline.hpp"
namespace geom {

template <class Point>
struct bezier_form : public bspline < Point >
{
    typedef bspline<Point> base_t;
    bezier_form(RAWTYPE(mk_stdvec(Point(0.0))) c):
        base_t(std::move(c),
                std::vector<double>(2 * c.size(), 0.0),
                c.size() - 1)
    {
        int deg =  base_t::degree();
        std::fill_n(base_t::t.begin() + (deg + 1), deg + 1, 1);
    }

    bezier_form(base_t && b):base_t(std::forward < base_t > (b))
    {}
};

template <class CptsT>
auto make_bezier_form(CptsT cpts) -> bezier_form<RAWTYPE(cpts[0])> 
{
	return bezier_form<RAWTYPE(cpts[0])>(std::move(cpts));
}

}

#endif // ASTI_BEZIER_FORM

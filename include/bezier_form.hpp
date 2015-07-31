#ifndef ASTI_BEZIER_FORM
#define ASTI_BEZIER_FORM
#include <vector>
#include "point_fwd.hpp"
namespace geom {

template <class Point>
struct bezier_form : public bspline < Point >
{

    bezier_form(RAWTYPE(mk_stdvec(Point(0.0))) c):
        bspline(std::move(c),
                std::vector<double>(2 * c.size(), 0.0),
                c.size() - 1)
    {
        int deg =  degree();
        std::fill_n(t.begin() + (deg + 1), deg + 1, 1);
    }

    bezier_form(bspline && b):bspline(std::forward < bspline > (b))
    {}
};

template <class CptsT>
auto make_bezier_form(CptsT cpts) -> bezier_form<RAWTYPE(cpts[0])> 
{
	return bezier_form<RAWTYPE(cpts[0])>(std::move(cpts));
}

}

#endif // ASTI_BEZIER_FORM

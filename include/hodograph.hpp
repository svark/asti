#ifndef ASTI_HODO_GRAPH
#define ASTI_HODO_GRAPH
#include "bspline_cons.hpp"
#include "bspline_queries.hpp"
#include <vector>
#include "type_utils.hpp"

namespace geom{

template <class Point>
bspline<Point> hodograph(const bspline<Point>& spl, int num_der = 1)
{
    std::vector<double> ks(spl.knots().begin() + num_der ,
                           spl.knots().end() - num_der );
    RAWTYPE(mk_stdvec(Point())) hodo_cpts(spl.control_points());
    auto const &cpts = spl.control_points();
    int d = spl.degree();
	auto const &t = spl.knots();

    for(int  k = 0 ; k < num_der; ++k ) {
		int i = k + 1;
        for(auto it = hodo_cpts.begin() + 1; it != hodo_cpts.end(); ++it, ++i)
        {
            auto prev_it = std::prev(it);
            *prev_it = make_pt((*it - *prev_it) * ((d-k)/t[i+d-k]-t[i]));
        }
    }
    hodo_cpts.erase(hodo_cpts.end() - num_der, hodo_cpts.end() );
    hodo_cpts.shrink_to_fit();
    return make_bspline(std::move(hodo_cpts),
                        std::move(ks) , d - num_der);
}
}
#endif // ASTI_HODO_GRAPH

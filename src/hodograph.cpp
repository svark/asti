#include "bspline_cons.hpp"
#include "hodograph.hpp"
#include "rmat.hpp"

namespace geom{

template <class Point>
bspline<Point> qry::hodograph(const bspline<Point>& spl,int num_der)
{
    std::vector<double> ks(spl.knots().begin() + num_der ,
                           spl.knots().end() - num_der );

    auto const &cpts = spl.control_points();
    int d = spl.degree();
	auto const &ts = spl.knots();

	ARRAY_TYPE(Point) hodo_cpts(cpts.size());

    rmat<Point>  mat(cpts, ts, d);
	auto hodo_cpts_iter = hodo_cpts.begin();

	const size_t inc = d - num_der + 1;
	for(size_t nu = d; nu < ts.size() - d - 1; nu += inc)
	{
		std::copy(cpts.begin() + nu - d, cpts.begin() + nu + 1, hodo_cpts_iter);
        mat.der_eval(nu, num_der, hodo_cpts_iter);
		hodo_cpts_iter += inc;
	}
	hodo_cpts.erase(hodo_cpts.end() - num_der, hodo_cpts.end());
    return make_bspline(std::move(hodo_cpts),
                        std::move(ks) , d - num_der);
}
}


#include "point.hpp"
/*
  Local Variables:
  eval:(load-file "./scripts/temp.el")
  eval:(instantiate-templates "hodograph" "qry" (list)  (list (cons "hodograph" (list  "double"
  "point2d_t" "point3d_t"  "point4d_t" ) )))
  End:
*/

namespace geom {
#include "hodograph_inst.inl"
}
//}}}

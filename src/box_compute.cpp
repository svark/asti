#include "point.hpp"
#include "box_compute.hpp"
#include "circle.hpp"
namespace geom {

template <class SplineType>
box<typename SplineType::point_t>
ops::compute_box(const SplineType &spl)
{
    box<typename SplineType::point_t> b;
    auto const & cpts = spl.control_points();
    for(auto const & c : cpts)
    {
        b += c;
    }
    return b;
}

template <class Point>
box<Point>
ops::compute_box(const circle<Point> &c)
{
    box<Point> b0;
    enum {dim = point_dim<Point>::dimension};
	typedef decltype(Point()-Point()) vector_t;
    for(int i =0; i < dim;++i)
    {
        Point p(0.0);
        p[i] = 1.0;

        vector_t ext_vec = cross (c.getPlaneNormal() ,
			decltype( c.getPlaneNormal()  )( make_vec(p)) );

        if(tol::eq(len(ext_vec), 0))
            continue;

        auto v = c.getRadius() * normalize(ext_vec);
        b0 += (c.getCenter() + v);
        b0 += (c.getCenter() - v);
    }
    return b0;
}

}
//{{{  instantiation scripts

/*
  Local Variables:
  eval:(load-file "./scripts/temp.el")
  eval:(setq methods (list "compute_box"
  ))
  eval:(setq spltypes (list "bspline<double>"
  "bspline<point2d_t>"
  "bspline<point3d_t>"
  "bspline<point4d_t>"
  "periodic_bspline<double>"
  "periodic_bspline<point2d_t>"
  "periodic_bspline<point3d_t>"
  "periodic_bspline<point4d_t>"
  "circle<point2d_t>"
  "circle<point3d_t>"
  ))
  eval:(instantiate-templates "box_compute" "ops" (list ) (list
  (cons (car methods) spltypes )))
  End:
// dump all explicitly instantiated templates below
*/
//}}}

//{{{  instantiation
#include "bspline.hpp"
#include "periodic_bspline.hpp"

namespace geom {
#include "box_compute_inst.inl"
}
//}}}

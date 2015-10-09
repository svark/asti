#include "box_compute.hpp"
#include "circle.hpp"
#include "conic.hpp"
#include "spline_traits.hpp"
#include "point_dim.hpp"
#include "bspline_x_cons.hpp"
#include "Eigen/Core"
#include "Eigen/LU"
#include "insert_knot.hpp"
#include "bspline_queries.hpp"
#include "smat.hpp"
#include <algorithm>
namespace geom {

template <class SplineType>
box<typename SplineType::point_t>
ops::compute_box(const SplineType &spl) {
    box<typename SplineType::point_t> b;
    auto const &cpts = spl.control_points();
    for (auto const &c : cpts) {
        b += c;
    }
    return b;
}

struct mintag : public std::false_type{};
struct maxtag : public std::true_type{};

namespace impl {

template <class SplineType, class MaxTag>
double
find_bound(MaxTag tag, SplineType &bs) {
    typedef typename SplineType::point_t point_t;
    static const int dimension = point_dim<point_t>::dimension;
    static_assert(dimension == 1, "1d spline is expected");
    const int d = bs.degree();
    while (true) {
        auto const &pts = bs.control_points();
        auto const &t = bs.knots();
        auto maxe = [&pts](){ return std::max_element(pts.begin(), pts.end());};
        auto mine = [&pts](){ return std::min_element(pts.begin(), pts.end());};
        auto pos = util::tag_switcher::eval(maxe, mine, tag);
        size_t i = std::distance(pts.begin(), pos);

        if (*pos == pts.front() || *pos == pts.back()) {
            return *pos;
        }
        auto pose = std::find_if(pos, pts.end(),
                                 [pos](double p){
                                     return tol::neqres(p,*pos);});

        if (pose == pts.end() || std::distance(pos,pose) >= d)
            return *pos;

        assert(!tol::eqres(*pos, *pose));
        // insert a knot u for which we would have pts[i] ==
        // pts[i+j]
        int j = pose - pos;

        Eigen::Matrix2d _2dx;
        _2dx << pts[i - 1] - pts[i], pts[i + j] - pts[i],
            t[i + d] - t[i]    , t[i + j] - t[i + d + j];

        Eigen::Vector2d _2dv;
        _2dv << pts[i + j] - pts[i], t[i + d] - t[i + d + j];
        _2dv = _2dx.lu().solve(_2dv);
		
        auto u = t[i + d] - _2dv[0] * (t[i + d] - t[i]);
        ops::insert_knot(bs, u).swap(bs);
    }
}
}

template <class SplineType>
box<typename SplineType::point_t>
ops::compute_box_tight(const SplineType &spl) {
    typedef typename SplineType::point_t point_t;
    typedef typename spline_traits<SplineType>::ptag ptag;
    typedef typename spline_traits<SplineType>::rtag rtag;
    const int dimension = point_dim<point_t>::dimension;

    box<point_t> b;
    int d = spl.degree();

    auto const &cpts = spl.control_points();
    for (int dim = 0; dim < dimension; ++dim) {
        std::vector<double> dimCpts;
        dimCpts.reserve(cpts.size());
        for (auto const &cp : cpts) {
            dimCpts.push_back(coord(cp, dim));
        }

        auto s = make_bspline(std::move(dimCpts), std::vector<double>(spl.knots()), d,
                              ptag(), rtag());
        auto bs( clamp_end(clamp_start(s)));

        coord_nonconst(b.lo, dim) = impl::find_bound(mintag(), bs);
        coord_nonconst(b.hi, dim) = impl::find_bound(maxtag(), bs);
    }

    return b;
}

template <class Point>
box<Point>
ops::compute_box(const circle<Point> &c) {
    box<Point> b0;
    enum { dim = point_dim<Point>::dimension };
    typedef decltype(Point() - Point()) vector_t;
    for (int i = 0; i < dim; ++i) {
        Point p(0.0);
        p[i] = 1.0;

        typedef decltype(c.plane_normal()) normal_t;
        vector_t ext_vec = cross(c.plane_normal(), normal_t(make_vec(p)));

        if (tol::small(len(ext_vec)))
            continue;

        auto v = c.radius() * normalize(ext_vec);
        b0 += (c.center() + v);
        b0 += (c.center() - v);
    }
    return b0;
}

template <class Point>
box<Point>
ops::compute_box(const conic_arc<Point> &c) {
    return compute_box(make_rbspline_from_conic(c));
}
}
//{{{  instantiation scripts

/*
  Local Variables:
  eval:(load-file "./scripts/temp.el")
  eval:(setq methods (list "compute_box"
  "compute_box_tight"
  ))
  eval:(setq spltypes (list "bspline<double>"
  "bspline<point2d_t>"
  "bspline<point3d_t>"
  "bspline<point4d_t>"
  "periodic_bspline<double>"
  "periodic_bspline<point2d_t>"
  "periodic_bspline<point3d_t>"
  "periodic_bspline<point4d_t>"
  ))
  eval:(setq alltypes  (apply 'append (list spltypes  (list "circle<point2d_t>"
  "circle<point3d_t>"
  "conic_arc<point3d_t>"
  "conic_arc<point2d_t>"))))
  eval:(instantiate-templates "box_compute" "ops" (list ) (list
  (cons (car methods) alltypes )
  (cons (cadr methods) spltypes )
  ))
  End:
  // dump all explicitly instantiated templates below
  */
//}}}

//{{{  instantiation

#include "point.hpp"
#include "periodic_bspline.hpp"

namespace geom {
#include "box_compute_inst.inl"
}
//}}}

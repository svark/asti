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
#include "tol.hpp"
namespace geom {

namespace impl{
template <class SplineType>
box<typename SplineType::point_t>
compute_box(const SplineType &spl, polynomial_tag);

template <class SplineType>
box<typename SplineType::point_t>
compute_box(const SplineType &spl, rational_tag);
}


template <class SplineType>
box<typename SplineType::point_t>
ops::compute_box(const SplineType &spl) {
    return impl::compute_box(spl, spline_traits<SplineType>::rtag());
}

template <class SplineType>
box<typename SplineType::point_t>
impl::compute_box(const SplineType &spl, polynomial_tag) {
    typedef typename SplineType::point_t point_t;
    box<point_t> b;
    auto const &cpts = spl.control_points();
    for (auto const &c : cpts) {
        b += c;
    }
    return b;
}

template <class SplineType>
box<typename SplineType::point_t>
impl::compute_box(const SplineType &spl, rational_tag) {
    typedef typename SplineType::point_t point_t;
    static const int dim = point_dim<point_t>::dimension;
    box<point_t> b;

    auto const &cpts = spl.control_points();
    typedef RAWTYPE(cpts[0]) pointw_t;

    auto comp = [](const pointw_t &p, const pointw_t& q){return  p[dim] < q[dim];};

    auto pos = std::max_element(cpts.begin(), cpts.end(), comp);
    assert(pos!=cpts.end());

    auto maxw = (*pos)[dim];
    assert(tol::not_small(maxw));

    for (auto const &c : cpts) {
        b += scaled_copy(lower_dim(c) , 1.0/maxw);
    }
    return b;
}

struct mintag : public std::false_type{};
struct maxtag : public std::true_type{};

namespace impl {
template <class SplineType, class MinMaxTag>
std::pair<double, bool>
find_knot_at_bound(const SplineType& bs, MinMaxTag tag)
{
    auto const &pts = bs.control_points();
    auto const &t   = bs.knots();
    const int d     = bs.degree();
    typedef RAWTYPE(pts[0]) pointw_t;

    auto comp = [](const pointw_t&p, const pointw_t &q)
    {
        return coord(p,0) < coord(q,0) - tol::resabs;
    };

    auto maxe = [&pts,&comp]() {
        return std::max_element(pts.begin(), pts.end(), comp);
    };

    auto mine = [&pts,&comp]() {
        return std::min_element(pts.begin(), pts.end(), comp);
    };

    auto pos = util::tag_switcher::eval(maxe, mine, tag);
    assert(pos != pts.end());

    size_t i = pos - pts.begin();
    double v = qry::greville(qry::get_spline(bs), i);

    if(*pos == pts.front() || *pos == pts.back()) {
        return std::make_pair(v, true);
    }

    auto neqres = [pos](const pointw_t& p) {
        return tol::neqres(coord(p,0), coord(*pos,0));
    };

    auto pose = std::find_if(pos, pts.end(), neqres);

    if (pose == pts.end() || pose - pos > d)
    {
        double vdash = qry::greville(qry::get_spline(bs), i + (d+1)/2);
        return std::make_pair(vdash, true);
    }

    assert(!tol::pt_eq(*pos, *pose));
    size_t j = pose - pos;


    // insert a knot u for which we would have pts[i] ==
    // pts[i+j] after insertion
    Eigen::Matrix2d _2dx;
    _2dx << coord(pts[i - 1],0) - coord(pts[i],0),
        coord(pts[i + j],0) - coord(pts[i],0),
            t[i + d]   - t[i]  , t[i + j]   - t[i + d + j];

    Eigen::FullPivLU<Eigen::Matrix2d> lu(_2dx);
    assert(lu.rank() == 2);

    Eigen::Vector2d _2dv;
    _2dv << coord(pts[i + j],0) - coord(pts[i],0),
        t[i + d] - t[i + d + j];
    _2dv = lu.solve(_2dv);

    auto u = t[i + d] - _2dv[0] * (t[i + d] - t[i]);
    return std::make_pair(u, false);
}

template <class SplineType, class MaxTag>
double
find_bound_by_insertion(SplineType &bs, MaxTag tag) {
    typedef typename SplineType::point_t point_t;
    static const int dimension = point_dim<point_t>::dimension;
    bool found_bound = false;
    double u = bs.param_range().second;
    std::tie(u, found_bound) = find_knot_at_bound(bs, tag);
    while(!found_bound) {
        ops::insert_knot(bs, u).swap(bs);
        std::tie(u, found_bound) = find_knot_at_bound(bs, tag);
    }
    return ncoord(bs.eval(u),0);
}

template <class SplineType>
box<typename SplineType::point_t>
compute_box_tight(const SplineType &spl, polynomial_tag) {
    typedef typename SplineType::point_t point_t;
    typedef typename spline_traits<SplineType>::ptag ptag;
    typedef typename spline_traits<SplineType>::rtag rtag;
    const int dimension = point_dim<point_t>::dimension;

    box<point_t> b;
    int d = spl.degree();

    auto const &cpts = spl.control_points();
#pragma loop(hint_parallel(4))
    for (int dim = 0; dim < dimension; ++dim) {
        std::vector<double> dimCpts;
        dimCpts.reserve(cpts.size());

        for (auto const &cp : cpts) {
            dimCpts.push_back(coord(cp, dim));
        }

        auto s = make_bspline(std::move(dimCpts),
                              std::vector<double>(spl.knots()), d,
                              ptag(), rtag());

        auto bs(clamp_end(clamp_start(std::move(s))));

        coord_nonconst(b.lo, dim) = find_bound_by_insertion(bs, mintag());
        coord_nonconst(b.hi, dim) = find_bound_by_insertion(bs, maxtag());
    }
    return b;
}

template <class SplineType>
box<typename SplineType::point_t>
compute_box_tight(const SplineType &spl, rational_tag) {
    typedef typename SplineType::point_t point_t;
    typedef typename spline_traits<SplineType>::ptag ptag;
    typedef typename spline_traits<SplineType>::rtag rtag;
    const int dimension = point_dim<point_t>::dimension;

    box<point_t> b;
    int d = spl.degree();

    auto const &cpts = spl.control_points();
#pragma loop(hint_parallel(4))
    for (int dim = 0; dim < dimension; ++dim) {
        ARRAY_TYPE(point2d_t) dimCpts;
        dimCpts.reserve(cpts.size());
        for (auto const &cp : cpts) {
            dimCpts.push_back(
                make_pt(coord(cp, dim),
                        coord(cp, dimension)));
        }

        auto s = make_bspline(std::move(dimCpts),
                              std::vector<double>(spl.knots()), d);

        auto bs(clamp_end(clamp_start(std::move(s))));

        coord_nonconst(b.lo, dim) = find_bound_by_insertion(bs, mintag());
        coord_nonconst(b.hi, dim) = find_bound_by_insertion(bs, maxtag());
    }

    return b;
}
}

template <class SplineType>
box<typename SplineType::point_t>
ops::compute_box_tight(const SplineType &spl) {
    typedef typename SplineType::point_t point_t;
    typedef typename spline_traits<SplineType>::rtag rtag;
    return impl::compute_box_tight(spl, rtag());
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
  "rational_bspline<double>"
  "rational_bspline<point2d_t>"
  "rational_bspline<point3d_t>"
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
#include "rational_bspline.hpp"
namespace geom {
#include "box_compute_inst.inl"
}
//}}}

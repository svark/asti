//-*- mode:c++ -*-
#ifndef ASTI_BSPLINE_SURFACE_H
#define ASTI_BSPLINE_SURFACE_H
#include <iterator>
#include "rmat.hpp"
#include "spline_traits.hpp"
#include "bspline.hpp"
#include <type_traits>
namespace geom {


enum bsurf_type {extruded_surf,ruled_surf, loft_surf, revolved_surf, generic};

template <class RTag = polynomial_tag,
          class PTagU = regular_tag,
          class PTagV = regular_tag,
          int surfType = generic >
struct bspline_surface_traits
{
    typedef RTag rtag;
    typedef PTagU ptagu;
    typedef PTagV ptagv;
    enum{ stype = surfType};
};

template <class Point, class SurfaceTraits =
          struct bspline_surface_traits<>  >
class bspline_surface {
public:

    typedef typename SurfaceTraits::ptagu ptagu;
    typedef typename SurfaceTraits::ptagv ptagv;
    typedef typename SurfaceTraits::rtag rtag;
    enum{is_rational = std::is_same<rtag, rational_tag>::value };

    typedef typename get_traits_type_from_tags<ptagu, rtag, Point>::type spl_traits;

    typedef Point point_t;
    typedef RAWTYPE(make_vec(point_t())) vec_t;

    typedef typename spl_traits::spline_type spl_t;

    typedef typename spl_t::point_t  pointw_t;
    typedef typename spl_t::vector_t vectorw_t;
    typedef typename spl_t::knots_t  knots_t;

    typedef typename spl_t::cpts_t  cpts_t;
    typedef typename spl_t::vcpts_t vwcpts_t;

    bspline_surface(cpts_t cpts_,
                    size_t stride_,
                    knots_t  t_u_,
                    knots_t  t_v_,
                    int degu_,
                    int degv_):stride(stride_),
                               cpts(std::move(cpts_)),
                               t_u(std::move(t_u_)),
                               t_v(std::move(t_v_)),
                               degu(degu_),
                               degv(degv_)
    {
    }

    point_t eval(double u, double v) const
    {
        return make_pt(eval_derivative(0, 0, u, v));
    }

    vec_t eval_derivative(int derOrderU,
                          int derOrderV,
                          double u, double v) const
    {
        u =  process_param(u, ptagu());
        v =  process_param(v, ptagv());

        assert(derOrderU >= 0 && derOrderV >= 0);
        return process_vec(u, v, derOrderU,
                           derOrderV, rtag());
    }

    std::pair<double, double> param_rangeu() const
    {
        return std::make_pair(t_u[degu], t_u[sizeu()]);
    }

    std::pair<double, double> param_rangev() const
    {
        return std::make_pair(t_v[degv], t_v[sizev()]);
    }

    const knots_t& knotsu() const { return t_u;}
    const knots_t& knotsv() const { return t_v;}

    const cpts_t& control_points() const { return cpts;}

private:
    double process_param(double u,  periodic_tag) const
    {
        return periodic_param(param_rangeu(u), u);
    }

    double process_param(double u,  regular_tag) const
    {
        return u;
    }

    vec_t process_vec(double u, double v,
                      int derOrderU,
                      int derOrderV,
                      rational_tag) const
    {
        vwcpts_t  dsu;
        for(int i = 0;i <= derOrderU; ++i)
        {
            vwcpts_t  dsv;
            for(int j = 0;j <= derOrderV; ++j)
            {
                dsv.push_back(process_vec(
                                  u, v, i,
                                  j, polynomial_tag()));
            }
            dsu.push_back(rational_derivatives(dsv).back());
        }
        return rational_derivatives(dsu).back();
    }

    vec_t process_vec(double u,
                      double v,
                      int    derOrderU,
                      int    derOrderV,
                      polynomial_tag) const
    {
        // as we are dealing with tensor product splines,  the
        // u and v derivatives can be computed independently
        // ./media/derivative-tpsurf.png

        std::vector<double> ub, vb;
        size_t nu_u, nu_v;

        std::tie(ub, nu_u) = rmat_base_vd(t_u, degu)
            .get_basis(u, derOrderU);

        std::tie(vb, nu_v) = rmat_base_vd(t_v, degv)
            .get_basis(v, derOrderV);

        vec_t base(0.0);

        int kj = 0;
        for(size_t j = nu_v - degv;j <= nu_v; ++j,++kj)
        {
            vec_t b(0.0);
            int ki = 0;
            for(size_t i = nu_u;i <= nu_u - degu; ++i,++ki)
                b += ub[ki] * make_vec(cpts[i + j * stride]);

            base += vb[kj] * b;
        }
        return base;
    }

    cpts_t cpts;
    size_t stride;
    knots_t t_u;
    knots_t t_v;
    int degu, degv;
};




}
#endif

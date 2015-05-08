#ifndef ASTI_BSPLINE_SURFACE_H
#define ASTI_BSPLINE_SURFACE_H
#include <iterator>

namespace geom {

template <class Point, class RTag = polynomial_tag,
          class PTagU = regular_tag,
          class PTagV = regular_tag>
struct bspline_surface_traits
{
    typedef RTag rtag;
    typedef PTagU ptagu;
    typedef PTagV ptagv;
};

template <class Point, class SurfaceTraits =
          struct bspline_surface_traits<Point>  >
class bspline_surface {
public:

    typedef typename SurfaceTraits::ptagu ptagu;
    typedef typename SurfaceTraits::ptagv ptagv;
    typedef typename SurfaceTraits::rtag rtag;
    typedef Point point_t;
    typedef decltype(mk_stdvec(point_t()))  cpts_t;
    typedef decltype(mk_stdvec(vector_t())) vcpts_t;

    bspline_surface(cpts_t cpts_,
                    size_t stride,
                    knots_u_t  t_u_,
                    knots_v_t  t_v_,
                    int degu_,
                    int degv_):stride(stride_),
                               cpts_t(std::move(cpts_)),
                               t_u(std::move(t_u_)),
                               t_v(std::move(t_v_)),
                               degu(degu_),
                               degv(degv_)
    {
    }

    Point eval(double u, double v) const
    {
        return make_pt(eval_derivative(0, 0, u, v));
    }

    vec_t eval_derivative(int derOrderU,
                          int derOrderV,
                          double u, double v) const
    {
        u =  process_param(u, ptagu());
        v =  process_param(v, ptagv());

        assert(order >= 0);
        return process_vec(u, v, derOrderU, derOrderV, rtag());
    }

    std::pair<double, double> param_rangeu() const
    {
        return std::make_pair(t_u[degu], t_u[sizeu()]);
    }

    std::pair<double, double> param_rangev() const
    {
        return std::make_pair(t_v[degv], t_v[sizev()]);
    }

    size_t sizev() const {return cpts_.size() / stride;}
    size_t sizeu() const {return stride;}

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
                      int derOrderU, in derOrderV,
                      rational_tag)
    {
        point_iter_traits < point_t * >::VectorContT  dsu;
        for(int i = 0;i <= derOrderU; ++i)
        {
            point_iter_traits < point_t * >::VectorContT  dsv;
            for(int j = 0;j <= derOrderV; ++j)
            {
                dsv.push_back(process_vec(
                                  u, v, i,
                                  j, polynomial_tag));
            }
            dsu.push_back(rational_derivatives(dsv).back());
        }
        return rational_derivatives(dsu).back();
    }

    vec_t process_vec(double u, double v, int derOrderU,
                      int derOrderV,
                      polynomial_tag)
    {
        // because we are dealing with tensor product splines,  the
        // u and v derivatives can be computed independently
        // ./media/derivative-tpsurf.png
        typedef rmat_base_vd::accumulator<Point> ac_t;
        ac_t  acu = rmat_base_vd(t_u, degu).
            get_accumulator < Point > (u, derOrderU);

        ac_t  acv = rmat_base_vd(t_v, degv).
            get_accumulator < Point > (v, derOrderV);

        cpts_t tmp;
        for(size_t j = 0 ;j < sizev(); ++j)
        {
            // swap out
            ac_t acu_tmp(std::move(acu));
            for(size_t i = 0;i < sizeu(); ++i)
                acv_tmp.prod(cpts_[i + j * stride]);

            acv.prod(acu_tmp.get());
            acu.swap(acu_tmp); // swap in
        }
        return acv.get();
    }

    cpts_t cpts;
    size_t stride;
    knots_u_t t_u;
    knots_v_t t_v;
    int degu, degv;
};


template <class VecVecT, class RTag,  class PTagU, class PTagV
bspline_surface < point_t, RTag, PTagU, PTagV >
make_bspline_surface(VecVecT cpts_,
                     knots_u_t t_u_,
                     knots_v_t t_v_,
                     int degu, int degv
    )
{
    typedef typename VecVecT::value_type::value_type point_t;
    typedef typename point_iter_traits <
        point_t * >::PointContainerT cpts_t;

    typedef decltype(cpts_.begin()) iter_t;
    return bspline_surface < point_t, RTag, PTagU, PTagV >
        ( std::accumulate(
            std::move_iterator<iter_t>(cpts_.begin()),
            std::move_iterator<iter_t>(cpts_.end()),
            cpts_t() ),
          cpts_.size(),
          std::move(t_u_),
          std::move(t_v_),
          degu, degv );
}

}
#endif

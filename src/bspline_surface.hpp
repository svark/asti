#ifndef BSPLINE_SURFACE_H
#define BSPLINE_SURFACE_H
#include <iterator>

namespace geom {

template <class Point>
class bspline_surface {
public:
    typedef typename point_traits < Point >::PointContainerT cpts_t;

    bspline_surface(cpts_t && cpts_, int stride,
                    const knots_u_t & t_u_,
                    const knots_v_t & t_v_,
                    int degu_, int degv_)
        :t_u(t_u_), t_v(t_v_), stride(stride_), cpts_t(cpts_),
         degu(degu_), degv(degv_)
    {
    }

    Point eval(double u, double v) const
    {
        return eval_derivative(0, 0, u, v);
    }
    Point eval_derivative(int numDerU, int numDerV,
                          double u, double v) const
    {
        // because we are dealing with tensor product splines,  the
        // u and v derivatives can be computed independently
        // ./media/derivative-tpsurf.png
        assert(numDer >= 0);
        typedef rmat_base_vd::accumulator<Point> ac_t;
        ac_t  acu = rmat_base_vd(t_u, degu).
            get_accumulator < Point > (u, numDerU);

        ac_t  acv = rmat_base_vd(t_v, degv).
            get_accumulator < Point > (v, numDerV);

        cpts_t tmp;
        for(int i = 0 ;i < sizeu(); ++i)
        {
            // swap out
            ac_t acv_tmp(std::move(acv));
            for(int j = 0;j < sizev(); ++j)
                acv_tmp.prod(cpts_[stride * i + j]);

            acu.prod(acv_tmp.get());
            acv.swap(acv_tmp); // swap in
        }
        return acu.get();
    }

    int sizeu() const {return cpts_.size() / stride;}
    int sizev() const {return stride;}
private:
    cpts_t cpts;
    int stride;
    knots_u_t t_u;
    knots_v_t t_v;
    int degu, degv;
};


template <class VecVecT>
make_bspline_surface(VecVecT cpts_,
                     const knots_u_t & t_u_,
                     const knots_v_t & t_v_,
                     int degu, int degv
    )
{
    typedef typename point_traits <
        typename VecVecT::value_type::value_type>::PointContainerT cpts_t;
     typedef decltype(cpts_.begin()) iter_t;
     return bspline_surface( std::accumulate(
                 std::move_iterator<iter_t>(cpts_.begin()),
                 std::move_iterator<iter_t>(cpts_.end()),
                 cpts_t() ), cpts_.size(), t_u_, t_v_, degu, degv);
}
}
#endif

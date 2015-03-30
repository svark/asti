#ifndef RATIONAL_BSPLINE_HPP
#define RATIONAL_BSPLINE_HPP
#include "point.hpp"
#include <vector>
#include <type_traits>

namespace geom
{
template <class SplineType>
struct rational_bspline
{
public:
    enum {dimension =  SplineType::dimension  - 1};
    typedef decltype(lower_dim(SplineType::point_t())) point_t;

    typedef typename SplineType::point_t pointw_t;
    typedef typename SplineType::vector_t vectorw_t;

    typedef decltype(SplineType::rebound_type<point_t>()) lowdim_spl_t;
    typedef typename lowdim_spl_t::vcpts_t vcpts_t;
    typedef SplineType         spl_t;
    typedef decltype(point_t() - point_t())  vector_t;

    typedef typename spl_t::cpts_t wcpts_t;
    typedef typename spl_t::cpts_t cpts_t;
    typedef typename spl_t::vcpts_t wvcpts_t;

    typedef std::vector<double> knots_t;

    rational_bspline(wcpts_t pts,
                     knots_t ks, int degree_);

    rational_bspline(const rational_bspline& other);

    rational_bspline(spl_t&& other)
        : spl(std::forward<spl_t>(other))
    {
    }

    template <class knot_iter>
    point_t blossom_eval(knot_iter f) const;

    // ./media/rational_bspline_der.png
    vcpts_t  eval_derivatives(int numDer, double u) const
    {
        vcpts_t vs; vs.reserve( numDer+1);
        std::vector<double>  ws; ws.reserve(numDer+1);
        std::vector<double> bbasis(numDer + 1, 0.0);
        auto const & vecs = spl.eval_derivatives(numDer, u);
        bbasis[0] =  1.0;
        for(int i = 0; i <= numDer; ++i) {
            // update the pascals triangle.
            for(int j = i - 1;j >= 1; --j)
                bbasis[j] = bbasis[j] + bbasis[j - 1];
            bbasis[i] = 1;
            auto vec =  vecs[i];
            auto vd(lower_dim(vec));
            auto w = vec[dimension];
            ws.push_back(w);
            for(int j = 1; j <= i; ++j)  {
                vector_t u(vs[i-j]);
                u  *= bbasis[j] * ws[j];
                vd -= u;
            }
            vd *= 1 / vecs[0][dimension];
            vs.push_back(vd);
        }
        return vs;
    }

    vector_t eval_derivative(int numDer, double u) const
    {
        return eval_derivatives(numDer, u).back();
    }

    point_t eval(double u) const;

    std::pair<double,double> param_range() const{
        return spl.param_range();
    }

    static point_t project(const pointw_t &pt)  {
       return _project(pt, std::integral_constant<int,dimension>());
    }

    template <int d>
    static point_t _project(const pointw_t &pt,
                            std::integral_constant<int,d>)  {
        point_t res;
        for(int i =0; i < dimension; ++i)
            res[i] = pt[i] / pt[dimension];
        return res;
    }
    static double _project(const pointw_t &pt,
                           std::integral_constant<int,1>) {
        return ( pt[0] / pt[1] );
    }

    const spl_t spline() const { return spl; }

    void swap(rational_bspline& other) {
        spl.swap(other.spl);
    }

    const knots_t & knots() const { return spl.knots(); }
    const wcpts_t &  control_points() const { return spl.control_points(); }
    int degree() const { return spl.degree(); };
    const typename spl_t::vector_t& base_point() const {
        return spl.base_point();
    }
    rational_bspline& translate(const vectorw_t& t) {
        spl.translate(t);
        return *this;
    }
private:
    spl_t  spl;
};

template <class SplineType>
rational_bspline<typename std::decay<SplineType>::type>
make_rbspline(SplineType&& spl)
{
    return rational_bspline<typename std::decay<SplineType>::type >( std::forward<SplineType>(spl) );
}


}
#endif // RATIONAL_BSPLINE_HPP

#ifndef ASTI_RATIONAL_BSPLINE_HPP
#define ASTI_RATIONAL_BSPLINE_HPP
#include "point.hpp"
#include "point_dim.hpp"
#include <vector>
#include "type_utils.hpp"
#include "spline_traits.hpp"
namespace geom
{

template <class VecsT>
auto rational_derivatives(const VecsT & vecs)
	->  decltype(mk_stdvec(lower_dim(vecs[0])))
{
    int derOrder = vecs.size() - 1;
    typedef decltype(lower_dim(vecs[0])) vec_t;
    enum {dim = point_dim < vec_t >::dimension};

    decltype(mk_stdvec(vec_t())) vs;
    vs.reserve(derOrder+1);

	std::vector<double>  ws;
	ws.reserve(derOrder+1);

	std::vector<double> bbasis(derOrder + 1, 0.0);
    bbasis[0] =  1.0;

	for(int i = 0; i <= derOrder; ++i) {
        // update the pascals triangle.
        for(int j = i - 1;j >= 1; --j)
            bbasis[j] = bbasis[j] + bbasis[j - 1];
        bbasis[i] = 1;
        auto vec =  vecs[i];
        auto vd(lower_dim(vec));
        auto w = vec[dim];
        ws.push_back(w);
        for(int j = 1; j <= i; ++j)  {
            vec_t u(vs[i-j]);
            u  *= bbasis[j] * ws[j];
            vd -= u;
        }
        vd *= 1 / vecs[0][dim];
        vs.push_back(vd);
    }
    return vs;
}

template <class Point, class PTag>
struct rational_bspline
{
public:
    // underlying spline type is one dimension higher
    typedef typename get_traits_type_from_tags < PTag, polynomial_tag,
                                                 typename inc_dimension < Point >::type
                                                 >::type spl_traits;

    typedef typename get_traits_type_from_tags<PTag, polynomial_tag, 
					       Point>::type ldim_spl_traits;

    typedef typename spl_traits::spline_type spl_t;
    enum {dimension =  spl_t::dimension  - 1};
    typedef Point point_t;

    typedef typename spl_t::point_t pointw_t;
    typedef typename spl_t::vector_t vectorw_t;

    typedef typename ldim_spl_traits::spline_type ldim_spl_t;
    typedef typename ldim_spl_t::vcpts_t vcpts_t;
    typedef typename ldim_spl_t::vector_t vector_t;


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
    vcpts_t  eval_derivatives(int derOrder, double u) const
    {
        return rational_derivatives(spl.eval_derivatives(derOrder, u));
    }

    vector_t eval_derivative(int derOrder, double u) const
    {
        return eval_derivatives(derOrder, u).back();
    }

    point_t eval(double u) const;

    std::pair<double,double> param_range() const{
        return spl.param_range();
    }

    static point_t project(const pointw_t &pt)  {
       return _project(pt, std::integral_constant<int,dimension>());
    }

    const spl_t spline() const { return spl; }

    void swap(rational_bspline& other) {
        spl.swap(other.spl);
    }

    const knots_t & knots() const {
        return spl.knots();
    }

    const wcpts_t & control_points() const {
        return spl.control_points();
    }

    int degree() const { return spl.degree(); };

    const typename spl_t::vector_t& base_point() const {
        return spl.base_point();
    }
    // store cpts relative to cg
    void optimize() { spl.optimize(); }

    // store cpts relative to origin
    void deoptimize() { spl.deoptimize(); }

    rational_bspline& translate(const vectorw_t& t) {
        spl.translate(t);
        return *this;
    }
private:

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

    spl_t  spl;
};

template <class WSplineType>
struct rbspline_from_spline
{
  typedef typename dec_dimension<typename WSplineType::point_t>::type  point_t;
  typedef rational_bspline<point_t,
                           typename spline_traits < WSplineType >::ptag
                           >  type;
};

template <class SplineType>
typename rbspline_from_spline < SplineType >::type
make_rbspline(SplineType&& spl)
{
  typedef typename rbspline_from_spline < SplineType >::type rspl_t;
  return rspl_t(std::forward<SplineType>(spl));
}


}
#endif // ASTI_RATIONAL_BSPLINE_HPP

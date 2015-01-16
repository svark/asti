#include "bspline.hpp"
#include <utility>
#include <vector>
namespace geom
{
  template <class Point>
  struct periodic_bspline
  {
    typedef Point   point_t;
    typedef decltype(Point() - Point())  vector_t;

    typedef std::vector<point_t> cpts_t;
    typedef std::vector<double> knots_t;
    typedef std::tuple<cpts_t&&,knots_t&&, int&&> tuple_t;

    static std::tuple<cpts_t&&,knots_t&&, int&&>
      wrap(const cpts_t& pts, const knots_t &ks, int degree);

    periodic_bspline(const cpts_t& pts, const knots_t &ks, int degree_)
      :spl(wrap(pts,ks,degree_))
    {
    }
    periodic_bspline(const periodic_bspline& other)
      :spl(other.spl)
    {
    }

    periodic_bspline(periodic_bspline&& other)
      :spl(std::forward<bspline<point_t>> (other.spl))
    {
    }

    double periodic_param(double u) const
    {
      double s,e;
      std::tie(s,e) = param_range();
      return s + fmod(u - s, (e-s));
    }

    point_t eval(double u) const {
      return spl.eval(periodic_param(u));
    }

    template <class knot_iter>
    point_t blossom_eval(knot_iter f) const {
      std::vector<double> modfs(degree);
      std::transform(f,
                     f + degree(),
                     modfs.begin(),
                     periodic_param);
      return spl.blossom_eval(modfs.begin());
    }

    vector_t eval_derivative(int numDer,double u) const
    {
      return spl.eval_derivative(numDer,periodic_param(u));
    }

    std::pair<double,double> param_range() const {
      auto b = spl.knots().cbegin();
      auto e = spl.knots().cend() - 1;
      return std::make_pair(*(b + degree() - 1),
                            *(e - degree() + 1) );
    }

    periodic_bspline&
    translate(const vector_t& t) {
      spl.translate(t); return *this;
    }

    /*constexpr*/ bool is_periodic() const { return true; }

    //constexpr
    static int dimension()  {
      return point_traits<point_t*>::dim;
    }
    // store cpts relative to cg
    void optimize() { spl.optimize(); }

    void swap( periodic_bspline & other ) {
      spl.swap(other.spl);
    }

    // accessors
    const knots_t & knots() const { return spl.knots(); }
    const cpts_t &  control_points() const { return spl.control_points(); }
    int degree() const { return spl.degree(); };
    const vector_t& base_point() const { return spl.base_point(); }
  private:
    bspline<point_t> spl;
  };

 
  
}

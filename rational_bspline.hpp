namespace geom
{
  template <int dim>
  class rational_bspline
  {
    typedef bspline<dim+1>        spl_t;
    typedef pt_t<dim>             point_t;
    typedef vec_t<dim>            vector_t;
    typedef pt_t<dim+1>           pointw_t;
    typedef std::vector<pointw_t> cpts_t;
    typedef std::vector<double> knots_t;


    rational_bspline(const cpts_t& pts,
                     const knots_t &ks, int degree_);

    rational_bspline(cpts_t&& pts,
                     knots_t &&ks, int degree_);

    rational_bspline(const std::vector<point_t>& pts,
                     const std::vector<double>& weights,
                     const knots_t &ks, int degree_);

    rational_bspline(cpts_t&& pts,
                     knots_t &&ks, int degree_);


    rational_bspline(const rational_bspline& other);

    rational_bspline(rational_bspline&& other);

    rational_bspline(bspline<dim+1>&& other);

    template <class knot_iter>
    point_t blossom_eval(knot_iter f);

    vector_t eval_derivative(double u) const;

    point_t eval(double u) const;

    std::pair<double,double> param_range() { return spl.param_range(); }

    point_t project(const pointw_t &pt)
    {
      point_t res;
      for(int i =0; i < dim; ++i)
        res.p[i] = pt.p[i] / pt.p[dim];
      return res;
    }
  private:
    spl_t & spl;
    weights_t weights;
  };

  template <int dim>
  rational_bspline<dim>
  make_rbspline(bspline<dim+1>&& spl)
  {
    return rational_bspline<dim>( std::forward(spl) );
  }

  template <int dim>
  rational_bspline<dim>
  make_rbspline(const bspline<dim+1>& spl)
  {
    return rational_bspline<dim>(spl);
  }
}

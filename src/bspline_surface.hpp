#ifndef BSPLINE_SURFACE_H
#define BSPLINE_SURFACE_H
namespace geom {

    class bspline_surface {
    public:
        bspline_surface(MatrixXd cpts_,
                        const knots_u_t & t_u_,
                        const knots_v_t & t_v_)
            :t_u(t_u_), t_v(t_v_), cpts(cpts_){}
    private:
        cpts_t cpts;
        knots_u_t t_u;
        knots_v_t t_v;
    };

}
#endif

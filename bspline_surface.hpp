#ifndef BSPLINE_SURFACE_H
#define BSPLINE_SURFACE_H
namespace geom {

    class bspline_surface {
    public:
        bspline_surface(cpts_t cpts_,
                        knots_u_t& t_u_,
                        knots_v_t t_v_);
    private:
        cpts_t cpts;
        knots_u_t t_u;
        knots_v_t t_v;
    };

}
#endif

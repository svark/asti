#ifndef LINE_HPP
#define LINE_HPP

#include <algorithm>
namespace geom{

    template <dim>
    struct line
    {
        line(const pt_t<dim>& p1,
             const vec_t<dim>& v)
            :start(p1), dir(v)
        {
            assert(tol::eq(len(dir) , 1.0));
        }

        static line<dim>
        make_line(const pt_t<dim>& p1,
                  const pt_t<dim>& p2
            )
        {
            vec_t<dim> dir = unit(p2-p1);
            return line<dim>(p1, dir);
        }

        pt_t<dim> point(double u) const
        {
            return center + dir * u;
        }

    private:
        pt_t<dim>  start;
        vec_t<dim> dir;
    };
}

#endif

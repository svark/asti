#ifndef ASTI_BSPLINE_CONS_HPP
#define ASTI_BSPLINE_CONS_HPP

#include "bspline.hpp"
#include "type_utils.hpp"
namespace geom {
namespace impl {

template <class Point, class CptsT>
bspline<Point>
make_bspline(CptsT && pts, std::vector<double>&& ks, int degree_, 
             std::true_type/* type of CptsT == bspline::cpts_t */ )
{
    typedef RAWTYPE(pts[0])  point_t;
    static_assert(std::is_same<point_t,Point>::value," deduced and passed val should be same");
    typedef bspline < point_t > spl_t;
    typedef typename spl_t::cpts_t cpts_t;
    typedef typename spl_t::knots_t knots_t;

    return spl_t(
        std::forward < cpts_t > (pts),
        std::forward < std::vector<double> > (ks),
        degree_);
}

template <class Point,class CptsT>
bspline<Point>
make_bspline(CptsT && pts, std::vector<double>&& ks, int degree_, std::false_type )
{

    typedef RAWTYPE(pts[0])  point_t;
    static_assert(std::is_same<point_t,Point>::value," deduced and passed val should be same");
    typedef bspline < point_t > spl_t;
    typedef typename spl_t::cpts_t cpts_t;
    typedef typename spl_t::knots_t knots_t;
    return spl_t(
        std::move(cpts_t(pts)),
        std::forward < knots_t > (ks),
        degree_);
}
}

template <class CptsT>
auto
make_bspline( CptsT pts, std::vector<double> ks, int degree_) 
    -> bspline < RAWTYPE(pts[0]) >
{
    typedef RAWTYPE(pts[0])  point_t;
    typedef bspline < point_t > spl_t;
    typedef typename spl_t::cpts_t cpts_t;

    return impl::make_bspline<point_t>(
        std::move(pts),
        std::move(ks),
        degree_,
        std::is_same < cpts_t, CptsT >  ()
        );
}

template <class Point>
bspline <Point>
make_bspline_arr(const Point *pts, 
                 const Point *ptsEnd, 
                 const double * ks,
                 const double *ksEnd,
                 int degree_)
{
    typedef Point  point_t;
    typedef bspline < point_t > spl_t;
    typedef typename spl_t::cpts_t cpts_t;
    typedef typename spl_t::knots_t knots_t;

    return spl_t(
        std::move(cpts_t(pts,ptsEnd)),
        std::move(knots_t(ks,ksEnd)),
        degree_);
}


template <class CptsT, class KnotsT>
auto
make_bspline( CptsT pts, KnotsT ks, int degree_)
    -> bspline < RAWTYPE(pts[0]) >
{
    typedef RAWTYPE(pts[0])  point_t;
    typedef bspline < point_t > spl_t;
    typedef typename spl_t::cpts_t cpts_t;
    typedef typename spl_t::knots_t knots_t;

    return impl::make_bspline<point_t>(
        std::move(pts),
        std::move(ks),
        degree_,
        std::is_same < cpts_t, CptsT > ()
        );
}



template <class CptsT, class KnotsT>
auto
make_bspline(std::tuple<CptsT,KnotsT, int>&& initVal )
    -> bspline < RAWTYPE(std::get<0>(initVal)[0]) >
{
    return make_bspline( std::forward<CptsT>(std::get<0>(initVal)),
                         std::forward<KnotsT>(std::get<1>(initVal)),
                         std::get<2>(initVal) );
}

}
#endif // ASTI_BSPLINE_CONS_HPP

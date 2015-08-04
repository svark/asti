#ifndef ASTI_CIRCLE_HPP
#define ASTI_CIRCLE_HPP
#include "rational_bspline_cons.hpp"
#include "point_dim.hpp"
#include <type_traits>
#include "geom_exception.hpp"
#include <math.h>
#include "smat.hpp"
/*

 */
#include <algorithm>
namespace geom{
template <class Point>
struct circle
{
    enum {dim = point_dim<Point>::dimension};
    typedef Point point_t;
    typedef decltype(make_vec(point_t())) vector_t;


    template <class PointU>
    struct allow_if_dim_2
    {
        typename std::enable_if< dim == 2 && point_dim < PointU > ::dimension ==2,
                                 PointU const & >::type type;
    };

    template <class PointU>
    circle(const PointU& center_,
           const PointU& point_,
           RAWTYPE(make_vec(PointU())) const & ydir_)
        :center(center_), start_pt(point_), ydir(ydir_)
    {
        auto x = start_pt - center;
        assert(!tol::eq(len(x),0));
        ydir -= x*dot((start_pt-center),ydir)/len(x);
        ydir = normalize(ydir) * len(x);
    }


    template <class PointU>
    circle(const PointU& center_, const PointU & point_, 
           typename std::enable_if<point_dim < PointU >::dimension == 2, int>::type = 0)
        :center(center_), start_pt(point_)
    {
        auto x = start_pt - center;
        ydir = normalize( make_vec( - x[1], x[0]) ) * len(x);
    }

    double
    static foot_param(const circle &c,
                      const point_t& p)
    {
        auto x = normalize(c.start_pt - c.center);
        auto y = normalize(c.ydir);

        auto v = (p - c.center);

        if( sqlen(v) < tol::sqresabs)
            throw geom_exception(point_at_axis_error);

        return atan2(dot(v,y) *y , dot(v,x) * x);
    }

    template <class PointIter,class ParamIter>
    void
    static foot_param(const circle &c,
                      PointIter ps,
                      PointIter end, ParamIter out)
    {
        auto x = normalize(c.start_pt - c.center);
        auto y = normalize(c.ydir);

        for( ;ps!=end; ++ps, ++out ) {

            auto p =*ps;
            auto v = (p - c.center);

            if( sqlen(v) < tol::sqresabs)
                throw geom_exception(point_at_axis_error);

            *out =  atan2(dot(v,y) *y , dot(v,x) * x);
        }
    }
    //  (@file :file-name "media/circle2.png" :to "./media/circle2.png" :display "eval at param")
    point_t eval(double u) const
    {
        auto x = (start_pt - center);
        auto y = ydir;
        return center + x * cos(u) + y * sin(u);
    }

    template <class ParamIter ,class PointIter>
    void eval(ParamIter us, ParamIter end, PointIter out) const
    {
        auto x = (start_pt - center);
        auto y = ydir;
        for( ;us!=end;++us,++out) {
            auto u = *us;
            *out = center + x * cos(u) + y * sin(u);
        }
    }

    vector_t tangent(double u) const
    {
        auto x = (start_pt - center);
        auto y = ydir;
        return  - x * sin(u) + y * cos(u);
    }

    vector_t normal(double u) const
    {
        return eval(u) - center;
    }
    // accessors
    point_t  getCenter() const
    {
        return center;
    }
    point_t getStart() const
    {
        return start_pt;
    }
    
	decltype(cross(vector_t(),vector_t()))
    getPlaneNormal() const
    {
        auto const &x = (start_pt - center);
        auto const &y = ydir;
        return  normalize( cross(x,y) );
    }

    double getRadius() const
    {
        return len(start_pt - center);
    }
private:
    point_t  center;
    point_t  start_pt;
    vector_t ydir;
};


template <class Point>
static circle<Point>
make_circle(const Point& p1,
            const Point& p2,
            const Point& p3)
{
    typedef Point point_t;
    auto v1 = p1 - p2;
    const double eps = tol::resabs;

    double lv1 = sqlen(v1);
    if(lv1 < eps*eps)
        throw geom_exception(circle_too_small);

    auto v2 = p2 - p3;
    double lv2 = sqlen(v2);
    if(lv2 < eps*eps)
        throw geom_exception(circle_too_small);

    auto v3 = p3 - p1;
    double lv3 = sqlen(v3);

    double d = sqlen( cross( v1, v2) );
    auto center = dot(v1, -v3) / d*d;
    double lens[] = {lv1, lv2, lv3};
    decltype(v1) vs[] = {v1,v2,v3};

    size_t idx = std::max_element(lens,lens+3) - lens;

    v1 = vs[(idx + 1)%3];   lv1 = lens[(idx + 1)%3];
    v2 = vs[(idx + 2)%3];   lv2 = lens[(idx + 2)%3];
    v3 = vs[idx];           lv3 = lens[idx];

    point_t pts[] = {p1,p2,p3};
    std::rotate(pts, pts + (idx + 1)%3, pts + 3);
    auto normal  = cross(v1,v2);
    double denom = sqlen(normal);
    if(denom < eps*eps )
        throw geom_exception(circle_too_small);

    double alpha = lv2 * dot(v1,-v3) *0.5 / ( denom );
    double beta  = lv3 * dot(-v1,v2) *0.5 / ( denom );

    auto cp = lerp(alpha,beta,pts[0],pts[1],pts[2]);
	decltype(normal) xdir(pts[0] - cp);
	auto nnormal (normalize(normal));
    return circle<point_t>(cp, pts[0], decltype(v1)(cross(nnormal, xdir)));
}

template <class Point>
static rational_bspline< Point, regular_tag >
to_rational(const circle<Point>& circ)
{
    auto start_pt =  circ.getStart();
    auto center =  circ.getCenter();
    auto const & x = (start_pt - center);
    decltype(x)  y = cross(circ.getPlaneNormal(), decltype(circ.getPlaneNormal())(x) );

    double radius =  circ.getRadius();
	
    auto a = start_pt + 2 * y * radius  * cos( M_PI/6.0);
    auto c = start_pt - 2 * y * radius  * cos( M_PI/6.0);
    auto b = center -  2 * x * radius ;
    
	auto p = lerp(0.5,a,c);
    auto q = lerp(0.5,a,b);
    auto r = lerp(0.5,c,b);

    double weights[] = {1,0.5,1,0.5,1,0.5,1};
    Point cpts[] = {p, a   ,q    ,b    ,r   ,c, p};
    double ts[]  = {0, 0, 0, 2*M_PI/3,2*M_PI/3,4*M_PI/3,4*M_PI/3, 2*M_PI,2*M_PI,2*M_PI};
	auto spl = make_bspline( 
     	interleave( mk_stdvec(cpts,cpts+ sizeof(cpts)/sizeof(Point)),
                    std::vector<double>( weights , weights + sizeof(weights)/sizeof(double) ) ), 
		std::vector<double>( ts , ts + sizeof(ts)/sizeof(double) ),
		2 );
	double st[] = {-2*M_PI/3,0,0};
	double es[] = {2*M_PI, 2*M_PI, 8.0*M_PI/3 };
	rebase_at_start(spl, st ).swap(spl);
	rebase_at_end(spl,es).swap(spl);

	return make_rbspline(std::move(spl));
		
}


}




#endif

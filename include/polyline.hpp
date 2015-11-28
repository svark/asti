#ifndef ASTI_POLYLINE
#define ASTI_POLYLINE
#include "bspline.hpp"
#include <vector>
namespace geom
{

template <class Point>
class polyline : public bspline<Point>
{

    typedef typename bspline<Point>::cpts_t cpts_t;
    typedef typename bspline<Point>::knots_t knots_t;
    
    static
    std::vector<double> compute_knots(const cpts_t &pts)
    {
        size_t num_pts = pts.size();
        std::vector<double> ks(num_pts+2);
        ks[0] = 0;
        ks[num_pts + 1] = 1.0;
        int prev  = 0;
        for(int i = 0; i < num_pts; ++i)
        {
            ks[i+1] = len(pts[i]-pts[prev]);
            prev = i;
        }
        return ks;
    }

    static
    std::vector<double> compute_knots(size_t num_pts, double start, double end)
    {
        assert(num_pts > 0);
        std::vector<double> ks(num_pts+2);
        ks[0] = start;
        ks[num_pts + 1] = end;

        for(size_t i = 0; i < num_pts; ++i)
        {
            ks[i+1] = start + i*(end - start)/(num_pts-1);
        }
        return ks;
    }

    static auto make_tuple(cpts_t pts, double start, double end)
    {
        int sz = pts.size();
        auto && ks = compute_knots(sz, start, end);
        return std::make_tuple(std::move(pts), std::forward<knots_t>(ks),1);
    }

    static auto make_tuple(cpts_t pts)
    {
        auto && ks = compute_knots(pts);
        return std::make_tuple(std::move(pts), std::forward<knots_t>(ks),1);
    }

public:
    polyline(cpts_t pts):
        bspline<Point>(make_tuple(std::move(cpts)))
    {
    }

    polyline(cpts_t pts, double start, double end):
        bspline<Point>(make_tuple(std::move(pts), start, end))
    {
    }
};

}
#endif // ASTI_POLYLINE

#ifndef ASTI_POLYLINE
#define ASTI_POLYLINE
#include "bspline.hpp"
#include <vector>
namespace geom
{

template <class Point>
class polyline : public bspline<Point>
{
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
        int prev  = 0;
        for(size_t i = 0; i < num_pts; ++i)
        {
            ks[i+1] = start + i*(end - start)/(num_pts-1);
            prev = i;
        }
        return ks;
    }
public:
    polyline(cpts_t pts):
        bspline<Point>(std::move(pts),compute_knots(pts),1)
    {
    }

    polyline(cpts_t pts, double start, double end):
        bspline<Point>(std::move(pts),compute_knots(pts.size(),start,end),1)
    {
    }
};

}
#endif // ASTI_POLYLINE

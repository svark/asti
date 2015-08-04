//-*- mode:c++ -*- 
#ifndef ASTI_TRIM_EXTEND_JOIN
#define ASTI_TRIM_EXTEND_JOIN
namespace geom
{
namespace ops
{
template <class SplineType>
extern SplineType
trim_curve(const SplineType &spl, double a, double b);

template <class SplineType>
extern SplineType
extract_regular_curve(const SplineType &spl);

template <class SplineType>
extern SplineType
extend_curve_start(const SplineType & spl, double delta);

template <class SplineType>
extern SplineType
extend_curve_end(const SplineType & spl, double delta);

template <class SplineType>
extern SplineType
extend_curve_end_to_pt(const SplineType & spl,
                       typename SplineType::point_t const & target);

template <class SplineType>
extern SplineType
join_starts(const SplineType& spl1,
            const SplineType& spl2,
            int join_cont);
}
}
#endif // ASTI_TRIM_EXTEND_JOIN

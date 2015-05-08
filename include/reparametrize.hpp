#pragma once
namespace geom {
namespace ops {

template <class SplineType>
extern SplineType reparametrize(const SplineType& spl,
                                double t1 = 0, double t2 = 1);
template <class SplineType>
extern SplineType reparametrize_start(const SplineType& spl,
                                      double t1 = 0);



}
}


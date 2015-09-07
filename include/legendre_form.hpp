#ifndef ASTI_LEGENDRE_FORM
#define ASTI_LEGENDRE_FORM
namespace geom {

template <class Point>
struct legendre_form
{
    typedef RAWTYPE(mk_stdvec(Point(0.0))) cpts_t;

    legendre_form(cpts_t a_,
                  double s_ = 0,
                  double e_ = 1.0)
        :a(std::move(a_)), s(s_), e(e_){}

    Point eval(double u) const{
        assert(a.size());
        Point v(0.0);
        double p0 = 1;
        u =  (u - s) / (e - s);
        u = 2 * u - 1;
        double p1 = u;
        v += p0 *a[0];
        v += p1 *a[1]*sqrt(3.0);
        for(size_t i = 2;i < a.size(); ++i)
        {
            auto pnu = ( (2 * i - 1) * u *  p1 -  (i - 1)*  p0 ) / i;
            v += pnu * a[i] * sqrt(2*i+1.0);
            p0 =  p1;
            p1 =  pnu;
        }
        return v;
    }

    Point  operator[](int i) const { return a[i];}
    size_t size() const            { return a.size(); }
    int    degree() const          { return int(a.size() - 1);}
    cpts_t const & coeffs() const  { return a;}
    Point  coeffs(int i) const  { return a[i];}
    double start_param() const     { return s;}
    double end_param() const       { return e;}
private:
    cpts_t a;
    double s, e;
};

}
#endif

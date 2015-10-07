#include "point.hpp"
#include <vector>
#include "catch.hpp"
#include "testutils.hpp"
#include "change_basis.hpp"
#include <algorithm>
#include "monomial_form.hpp"
#include "bezier_form.hpp"
#include "legendre_form.hpp"

using geom::monomial_form;
using geom::bezier_form;
using geom::ops::to_bezier;
using geom::ops::to_monomial;
using geom::ops::to_legendre;
using geom::legendre_form;
TEST_CASE("change basis", "[bezier][monomial]"){
    SECTION("monomial to bezier") {
        std::vector<double> mon(5);
        for(int i = 0;i < 5; ++i)
            mon[i] = i + 1;

        monomial_form < double >  mf(mon, 0, 1.0);
        auto const &bzf      = to_bezier(mf);
        REQUIRE(bzf.eval(.1) == Approx(mf.eval(.1)));
        auto const &mf_dual  = to_monomial(bzf);
        auto const &cfs      = mf_dual.coeffs();
        REQUIRE(cfs.size()  == mon.size());
        for(size_t i = 0;i < cfs.size(); ++i)
        {
            REQUIRE(cfs[i] == Approx(mon[i]));
        }
    }
    SECTION("bezier  to monomial") {
        std::vector<double> c(5);
        for(int i = 0;i < 5; ++i)
            c[i] = i + 1;

        bezier_form < double > bf(c, 1, 2.0);
        auto const & bf_dual = to_bezier (to_monomial(bf));
        auto const &  cfs = bf_dual.control_points();
        REQUIRE(cfs.size() == c.size());
        for(size_t i = 0;i < cfs.size(); ++i)
        {
            REQUIRE(cfs[i] == Approx(c[i]));
        }
    }
    SECTION("legendre to bezier") {
        std::vector<double> legf(5);
        for(int i = 0;i < 5; ++i)
            legf[i] = i + 1;

        legendre_form < double >  lf(legf, 0, 1.0);
        auto const &bzf      = to_bezier(lf);
        REQUIRE(bzf.eval(.1) == Approx(lf.eval(.1)));
        REQUIRE(bzf.eval(.9) == Approx(lf.eval(.9)));
        auto const &lf_dual  = to_legendre(bzf);
        auto const &cfs      = lf_dual.coeffs();
        REQUIRE(cfs.size()  == legf.size());
        for(size_t i = 0;i < cfs.size(); ++i)
        {
            REQUIRE(cfs[i] == Approx(legf[i]));
        }
    }
    SECTION("bezier to legendre") {
        std::vector<double> c(5);
        for(int i = 0;i < 5; ++i)
            c[i] = i + 1;

        bezier_form < double > bf(c, 1, 2.0);
        auto const & lf = to_legendre(bf);
        auto const & bf_dual = to_bezier(lf);
        REQUIRE(bf.eval(1.1) == Approx(lf.eval(1.1)));
        REQUIRE(bf.eval(1.9) == Approx(lf.eval(1.9)));
        auto const &  cfs = bf_dual.control_points();
        REQUIRE(cfs.size() == c.size());
        for(size_t i = 0;i < cfs.size(); ++i)
        {
            REQUIRE(cfs[i] == Approx(c[i]));
        }
    }
}

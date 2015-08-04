#include "point.hpp"
#include <vector>
#include "catch.hpp"
#include "testutils.hpp"
#include "change_basis.hpp"
#include <algorithm>
#include "monomial_form.hpp"
#include "bezier_form.hpp"

using geom::monomial_form;
using geom::bezier_form;
using geom::ops::to_bezier;
using geom::ops::to_monomial;
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
}

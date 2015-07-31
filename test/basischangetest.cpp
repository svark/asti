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

        monomial_form < double >  mf(mon);
		auto const &mf_dual = to_monomial(to_bezier(mf));
        auto const &cfs = mf_dual.coeffs();
        REQUIRE(cfs.size() == mon.size());
        for(size_t i = 0;i < cfs.size(); ++i)
        {
            REQUIRE(cfs[i] == Approx(mon[i]));
        }
    }
    SECTION("bezier  to monomial") {
        std::vector<double> c(5);
        for(int i = 0;i < 5; ++i)
            c[i] = i + 1;

        bezier_form < double > bf(c);
		auto & bf_dual = to_bezier (to_monomial(bf));
        auto &  cfs = bf_dual.control_points();
        REQUIRE(cfs.size() == c.size());
        for(size_t i = 0;i < cfs.size(); ++i)
        {
            REQUIRE(cfs[i] == Approx(c[i]));
        }
    }
}

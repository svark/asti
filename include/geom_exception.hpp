#ifndef ASTI_GEOM_EXCEPTION_HPP
#define ASTI_GEOM_EXCEPTION_HPP
#include <string>
#include <map>
namespace geom {
enum geom_error_code_t
{
    no_error = 0,
    knot_not_in_range_error = 1,
    knot_not_in_range_error_der = 2,
    bad_knot_spacing_t = 3,
    unknown_interp_option_t= 4,
    interp_ends_dont_match_for_periodic = 5,
    continuity_condition_too_tight= 6,
    vectors_not_in_plane_of_points = 7,
    lines_do_not_meet = 8,
    point_at_axis_error = 9,
    circle_too_small = 10,
    knots_incompatible_for_merging = 11,
    tangent_vectors_too_small = 12,
    degenerate_or_small_conic = 13,
    bspline_invariants_violated = 14,
    mismatched_array_sizes = 15,
    invalid_periodic_data = 16,
    duplicate_point_data = 17,
    degenerate_circle = 18,
	circle_too_large = 19
};

extern  std::map<geom_error_code_t, std::string> init_dict();

struct geom_exception
{
    geom_exception(geom_error_code_t code, const char *func, int l)
        :exception_code(code), function(func), line(l)
    {
    }

    geom_exception(geom_error_code_t code)
        :exception_code(code),line(-1)
    {
    }

    std::string what() { return dict().at(code());}
    geom_error_code_t code() { return exception_code;}
private:
    static const std::map<geom_error_code_t, std::string> & dict() {
        static std::map<geom_error_code_t, std::string> dict_ ( std::move(init_dict()));
        return dict_;
    }
    geom_error_code_t exception_code;
    std::string function;
    int line;
};
}
#endif

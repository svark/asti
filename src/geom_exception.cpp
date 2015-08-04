#include "geom_exception.hpp"

namespace geom {

std::map<geom_error_code_t, std::string>
init_dict()
{
    std::map<geom_error_code_t, std::string> dict_;
    auto add =  [&dict_](geom_error_code_t x, const char *y) {
        dict_.insert(std::make_pair(x, std::string(y)));
    };

    add(knot_not_in_range_error,"parameter is not in knot range");
    add(knot_not_in_range_error_der,"spline is not smooth enough at the given parameter for derivatives to be computed");
    add(bad_knot_spacing_t,
        "knot range is too close to reparametrise");

    add(unknown_interp_option_t, "cannot recognize parameterization option. \nSee documentation for options for parametrization");
    add(interp_ends_dont_match_for_periodic, "Point ends do not match for periodic interpolation");
    add(continuity_condition_too_tight,"Join continutiy must not exceed or equal degree");
    add(vectors_not_in_plane_of_points,"Vectors provided as tangents to conic are not the plane of given points");
    add(lines_do_not_meet,"Given lines do not meet");
    add(circle_too_small, "Circle is tool small compared to precision of computation");
    add(point_at_axis_error, "Nearest point computation failed.\nGiven point is very near the axis of circle");
    add(knots_incompatible_for_merging,"incompatible knots in given curves");
    return dict_;
}

}

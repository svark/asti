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
    add(point_at_axis_error, "Nearest point computation failed.\n The given point is very near the axis of circle");
    add(knots_incompatible_for_merging,"incompatible knots in given curves");
    add(tangent_vectors_too_small,"tangent at the start or end are too small in magnitude");

    add(degenerate_or_small_conic,"points defining the conic would form a degenerate or too small a conic");
    add(bspline_invariants_violated, "bspline invariants have been violated, this could one of the following"
        "\n 1) degree non negative\n"
        " 2) unique knots must be atleast 2\n"
        " 3)control points should be atleast as many as degree\n"
        " 4)num_cpts + degree + 1 == num_knots\n" );

    add( mismatched_array_sizes, "spline interpolation expectes that number of tangents provided matches number of points,\n"
         "note: tangents may be zero in which case they will be ignored" );
    add( invalid_periodic_data ,"end point must match start point for periodic spline interpolation");
    add( duplicate_point_data ,"points to interpolate may not have consecutive duplicates");
    add( degenerate_circle ,"circle generated would be degenerate or small");
    return dict_;
}

}

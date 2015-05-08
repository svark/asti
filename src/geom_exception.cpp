#include "geom_exception.hpp"

namespace geom {



std::map<geom_error_code_t, std::string>
init_dict()
  {
        std::map<geom_error_code_t, std::string> dict_;
      dict_.insert(std::make_pair(knot_not_in_range_error,
                                  std::string("parameter is not in knot range")));
      dict_.insert(std::make_pair(knot_not_in_range_error_der,
                                  std::string("spline is not smooth enough at the given parameter for derivatives to be computed")));
      dict_.insert(std::make_pair(bad_knot_spacing_t, std::string("knot range is too close to reparametrise")));

      dict_.insert(std::make_pair(unknown_interp_option_t, std::string("cannot recognize parameterization option. \nSee documentation for options for parametrization")));
      dict_.insert(std::make_pair(interp_ends_dont_match_for_periodic, std::string("Point ends do not match for periodic interpolation")));
      dict_.insert(std::make_pair(continuity_condition_too_tight,std::string("Join continutiy must not exceed or equal degree")));
      dict_.insert(std::make_pair(vectors_not_in_plane_of_points,std::string("Vectors provided as tangents to conic are not the plane of given points")));
      dict_.insert(std::make_pair(lines_do_not_meet,std::string("Given lines do not meet")));
      dict_.insert(std::make_pair(circle_too_small, std::string("Circle is tool small compared to precision of computation")));
      dict_.insert(std::make_pair(point_at_axis_error, std::string("Nearest point computation failed.\nGiven point is very near the axis of circle")));
      dict_.insert(std::make_pair(knots_incompatible_for_merging,
                                    std::string("incompatible knots in given curves")));
      return dict_;
   }

}

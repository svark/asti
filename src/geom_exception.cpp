#include "stdafx.h"
#include "geom_exception.hpp"

namespace geom {

	

std::map<geom_error_code_t, std::string> 
init_dict()
  {
	    std::map<geom_error_code_t, std::string> dict_;
      dict_.insert(std::make_pair(knot_not_in_range_error, std::string("parameter is not in knot range")));
      dict_.insert(std::make_pair(knot_not_in_range_error_der, std::string("spline is not smooth enough for derivatives")));
	  dict_.insert(std::make_pair(bad_knot_spacing_t, std::string("knot range is too close to reparametrise")));
	  return dict_;
   }

}
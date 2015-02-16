#include <map>
namespace geom {
   enum geom_error_code_t
   {
      knot_not_in_range_error = 1,
      knot_not_in_range_error_der = 2,
      spline_has_no_root = 3,
	  bad_knot_spacing_t = 4,
   };
   
   extern  std::map<geom_error_code_t, std::string> init_dict();

   struct geom_exception
   {
       geom_exception(geom_error_code_t code, const char *func, int l)
       :code_(code), function(func), line(l)
      {
      }

      geom_exception(geom_error_code_t code)
        :code_(code),line(-1)
      {
      }
	  
      std::string what() { return dict().at(code());}
      geom_error_code_t code() { return code_;}
   private:
	   static const std::map<geom_error_code_t, std::string> & dict() {
		  static std::map<geom_error_code_t, std::string> dict_ ( std::move(init_dict()));
		  return dict_;
	  }
      geom_error_code_t code_;
      std::string function;
      int line;
   };
}

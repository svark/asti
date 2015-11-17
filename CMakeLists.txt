if (NOT FIPS_IMPORT)
    cmake_minimum_required(VERSION 2.8)
    get_filename_component(FIPS_ROOT_DIR "../fips" ABSOLUTE)
    include("${FIPS_ROOT_DIR}/cmake/fips.cmake")
    set(FIPS_EXCEPTIONS ON CACHE BOOL "Enable C++ exceptions" FORCE)
    fips_setup()
    fips_project(asti)
    message (STATUS "fips setup complete")
endif()


if( LINUX )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -gnu=c++1y")
endif()
IF(APPLE)
   set(CPP_11_FLAGS "-std=c++11 -stdlib=libc++")
ELSEIF(UNIX)
   set(CPP_11_FLAGS "-std=c++11")
ENDIF(APPLE)

fips_begin_lib(asti)
   include_directories(include)

   file(GLOB asti_sources src/*.cpp)
   fips_files(${asti_sources})
fips_end_lib()

message (STATUS "fips asti lib build complete")

fips_begin_app(test_asti cmdline)
    include_directories(include)
    file(GLOB asti_test_files test/*.cpp)
    fips_files(${asti_test_files})
    fips_deps(asti)
    add_test(test_asti test_asti)
fips_end_app()


if (NOT FIPS_IMPORT)
    fips_finish()
endif()

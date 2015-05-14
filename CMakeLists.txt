# This CMakeLists.txt file helps defining your block building and compiling
# Include the main biicode macros and functions
# To learn more about the CMake use with biicode, visit http://docs.biicode.com/c++/building.html
# Or check the examples below


include(${CMAKE_HOME_DIRECTORY}/biicode.cmake)
message (STATUS "atttempting to include boost")
include(biicode/boost/setup)

# # Initializes block variables
#INIT_BIICODE_BLOCK()
# # This function creates the following variables:
# #     ${BII_BLOCK_NAME}       The name of the current block (e.g. "box2d")
# #     ${BII_BLOCK_USER}       The user's name (e.g. "phil")
# #     ${BII_BLOCK_PREFIX}     The directory where the block is located ("blocks" or "deps")

# # Also it loads variables from the cmake/bii_user_block_vars.cmake
# #     ${BII_CREATE_LIB}       TRUE if you want to create the library
# #     ${BII_LIB_SRC}          File list to create the library
# #     ${BII_LIB_TYPE}         STATIC(default) or SHARED
# #     ${BII_LIB_DATA_FILES}   Data files that have to be copied to bin
# #     ${BII_LIB_DEPS}         Dependencies to other libraries (user2_block2, user3_blockX)
# #     ${BII_LIB_SYSTEM_DEPS}  System linking requirements as winmm, m, ws32, pthread...

# # You can use or modify them here, for example, to add or remove files from targets based on OS
# # Or use typical cmake configurations done BEFORE defining targets. Examples:
# #     ADD_DEFINITIONS(-DFOO)
# #     FIND_PACKAGE(OpenGL QUIET)
# #     BII_FILTER_LIB_SRC(${BII_LIB_SRC})
# #     You can add INCLUDE_DIRECTORIES here too


#SET(BII_LIB_SRC ${SRC_LIST})
# # Actually create targets: EXEcutables and libraries.


file(GLOB files "test/*.cpp")
message (STATUS "${files}")
foreach(file ${files})
  file(RELATIVE_PATH RELP ${CMAKE_CURRENT_SOURCE_DIR} ${file})
  LIST(REMOVE_ITEM BII_LIB_SRC ${RELP})
endforeach()

foreach(file ${files})
  file(RELATIVE_PATH RELP ${CMAKE_CURRENT_SOURCE_DIR} ${file})
  LIST(APPEND BII_test_main_SRC ${RELP})
endforeach()

SET(BII_BLOCK_TESTS test_main)

ADD_BII_TARGETS()
SET(BII_CREATE_LIB True)
bii_find_boost()

IF(APPLE)
   set(CPP_11_FLAGS "-std=c++11 -stdlib=libc++")
ELSEIF(UNIX)
   set(CPP_11_FLAGS "-std=c++11")
ENDIF(APPLE)

if(Boost_FOUND)
    message(STATUS "found boost")
    target_compile_options(${BII_BLOCK_TARGET} INTERFACE ${CPP_11_FLAGS})
endif()

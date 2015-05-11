# This CMakeLists.txt file helps defining your block building and compiling
# Include the main biicode macros and functions
# To learn more about the CMake use with biicode, visit http://docs.biicode.com/c++/building.html
# Or check the examples below


include(${CMAKE_HOME_DIRECTORY}/biicode.cmake)
message (STATUS "atttempting to including boost")
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
message (${BII_LIB_SRC})
message(${BII_TEST_SRC})
# # Actually create targets: EXEcutables and libraries.
ADD_BII_TARGETS()

SET(BII_CREATE_LIB True)
# # This function creates the following variables:
# #     ${BII_BLOCK_TARGETS} List of targets defined in this block
# #     ${BII_LIB_TARGET}  Target library name, usually in the form "user_block"
# #     ${BII_exe_name_TARGET}: Executable target (e.g. ${BII_main_TARGET}. You can also use
# #                            directly the name of the executable target (e.g. user_block_main)

# # Example, if your library links with pthread:
# #     TARGET_LINK_LIBRARIES(${BII_LIB_TARGET} pthread)
# #        or directly use the library target name
# #     TARGET_LINK_LIBRARIES(user_block pthread)
# #     NOTE:  This can be also done adding pthread to ${BII_LIB_DEPS}
# #            BEFORE calling ADD_BIICODE_TARGETS()
# #
# # Example: how to activate C++11. Uncomment one of this
# #     target_compile_options(${BII_LIB_TARGET} PUBLIC -std=c++11) # Win, or linux if not gcc problems
# #     target_compile_options(${BII_LIB_TARGET} PUBLIC -std=c++11  -stdlib=libc++) # MacOS
# #     target_compile_options(${BII_LIB_TARGET} PUBLIC -std=c++11 -Wl,--no-as-needed) # Ubuntu if gcc problems
# # with LLVM, you might need to link the target with "-lc++abi"
# #
# # Note: If you're in Ubuntu or Windows OS and you get an error with -std=c++11 option , try to use -std=c++0x instead.
# #
# # If you want your code to be portable and activate properly C++11, you have to handle it:
# # if(WIN32)
# #     target_compile_options(${BII_LIB_TARGET} ... windows
# # elseif(APPLE)
# #     target_compile_options(${BII_LIB_TARGET} ... MacOS
# # elseif(...)
# #     ...
# # endif()


bii_find_boost()
IF(APPLE)
   set(CPP_11_FLAGS "-std=c++11 -stdlib=libc++")
ELSEIF (WIN32 OR UNIX)
   set(CPP_11_FLAGS "-std=c++11")
ENDIF(APPLE)

if(Boost_FOUND)
    message(STATUS "found boost")
    target_compile_options(${BII_BLOCK_TARGET} INTERFACE ${CPP_11_FLAGS})
endif()

project(JADE++)
cmake_minimum_required(VERSION 2.8.8)
set(CMAKE_MODULE_PATH ${JADE++_SOURCE_DIR}/cmake/modules ${CMAKE_MODULE_PATH})
include(MacroOutOfSourceBuild)
macro_ensure_out_of_source_build("${PROJECT_NAME} requires an out of source build.")
add_subdirectory(src)
#EXEC_PROGRAM(executable [args ...])

# If doxygen exists, use the doc/CMakeLists.txt to generate code documentation.
# find_package( Doxygen )
# if( DOXYGEN_FOUND )
#   # This is the directory that contains all the source files we want to have
#   # documented.
#   set( DOXYGEN_SOURCE_DIR ${OnzaFDTD_SOURCE_DIR} )
#   # This is the directory that holds the doxygen doxyfile template (doxy.conf.in)
#   set( DOXYGEN_INPUT_DIR  ${OnzaFDTD_SOURCE_DIR}/doc )
#   # This is the directory where the generated code documentation goes.
#   set( DOXYGEN_OUTPUT_DIR ${OnzaFDTD_BINARY_DIR}/doc )
#   # It is possible to put adddditional CMake commands in sub-directory
#   # CMakeLists.txt, but the add_subdirectory command must be used to specify
#   # which sub-directories should be scanned.
#   #add_subdirectory( doc )
# else()
#   message( WARNING "Doxygen not found - code documentation won't be created" )
# endif()


# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
file(GLOB current_dir_src *.cc) 
get_filename_component(lib_name ${CMAKE_CURRENT_SOURCE_DIR} NAME)
add_library(${lib_name} ${current_dir_src})
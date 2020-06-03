# set type specific output defaults
include(GNUInstallDirs)

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")

find_package(MPI REQUIRED)
IF(NOT MPI_CXX_FOUND)
    MESSAGE(FATAL_ERROR "\n"
        "Could not find any suitable mpi library!\n"
        ${MPI_ADDITIONAL_ERROR_STRING}
        "\nPlease ensure that an mpi library is installed on your computer\n"
        "and set CMAKE_CXX_COMPILER to the appropriate mpi wrappers:\n"
        "    $ CXX=\".../mpicxx\" cmake <...>\n"
        "    $ cmake -DCMAKE_CXX_COMPILER=\".../mpicxx\" <...>\n"
        "Or with additional C and Fortran wrappers (recommended!):\n"
        "    $ CC=\".../mpicc\" CXX=\".../mpicxx\" F90=\".../mpif90\" cmake <...>\n"
        "    $ cmake -DCMAKE_C_COMPILER=\".../mpicc\"\\\n"
        "            -DCMAKE_CXX_COMPILER=\".../mpicxx\"\\\n"
        "            -DCMAKE_Fortran_COMPILER=\".../mpif90\"\\\n"
        "            <...>\n"
        )
endif()
set(LIB JADEXX)
set(HEADER_LIST "include/jade.h")
add_library(${LIB}
        jade.cc
        ${HEADER_LIST}
)
add_library(${NAMESPACE}${LIB} ALIAS ${LIB})
target_include_directories(${LIB} PUBLIC include)
target_link_libraries(${LIB} PUBLIC MPI::MPI_CXX)
target_compile_features(${LIB} PUBLIC cxx_std_11)

# IDEs should put the headers in a nice place
source_group(TREE "${PROJECT_SOURCE_DIR}/src/lib/include" PREFIX "Header Files" FILES ${HEADER_LIST})
## installation and export of targets
#
#install(TARGETS ${LIB} EXPORT ${TARGETS_EXPORT_NAME}
#    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
#    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT RuntimeLibraries
#    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT RuntimeLibraries
#    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
#)

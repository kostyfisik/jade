include(GNUInstallDirs)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/tests/${CMAKE_INSTALL_BINDIR}")

function(add_jadexx_test test_src libraries)
    get_filename_component(test_name ${test_src} NAME_WE)

    add_executable(${test_name} ${test_src})
    foreach(lib IN LISTS libraries)
        message("linking ${test_name} with ${lib}")
        target_link_libraries(${test_name} ${lib})
    endforeach()
    add_test(NAME ${test_name} COMMAND ${test_name})
endfunction()

foreach(jadexx_test
        test-jade-single-function.cc
        test-jade-feed.cc
        test-jade.cc
)
    set(target_libs JADEXX)
    add_jadexx_test(${jadexx_test} "${target_libs}")
endforeach()

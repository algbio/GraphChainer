# Download and unpack lemon at configure time
configure_file(CMakeLists.txt.in lemon-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lemon-download )
if(result)
    message(FATAL_ERROR "CMake step for lemon failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lemon-download )
if(result)
    message(FATAL_ERROR "Build step for lemon failed: ${result}")
endif()


# Add lemon directly to our build. This defines
# the lemon target.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/lemon-src
        ${CMAKE_CURRENT_BINARY_DIR}/lemon-build
        EXCLUDE_FROM_ALL)
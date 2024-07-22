# #[[

# @file   CMakeLists.txt
# @brief  a test executable 

# ]]

file(GLOB_RECURSE TEST_FILES "*.cpp")

# Add test executables
foreach(test_file ${TEST_FILES})
    get_filename_component(test_name ${test_file} NAME_WE)
    string(REPLACE "-test" "" test_name ${test_name}) 
    add_executable(${test_name} ${test_file})
    target_include_directories(${test_name} PUBLIC ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/tests)
    if(APPLE)
        target_link_libraries(${test_name} PRIVATE "/opt/homebrew/opt/libomp/lib/libomp.dylib" ${PROJECT_NAME})
    else()
        target_link_libraries(${test_name} PRIVATE OpenMP::OpenMP_CXX ${PROJECT_NAME})
    endif()
    target_compile_options(${test_name} PRIVATE -Wall -Wpedantic -O0)
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()
#[[
@file   CMakeLists.txt
@brief  a test executable 

]]

#################################
# make executables for test.cpp #
#################################

add_executable(SIMPLEPIR-test test-simplepir.cpp)
target_compile_options(SIMPLEPIR-test PRIVATE -std=c++11 -O3)
target_link_libraries(SIMPLEPIR-test PRIVATE simplepir)

add_executable(NTT-test test-ntt.cpp)
target_compile_options(NTT-test PRIVATE -std=c++11 -O3)
target_link_libraries(NTT-test PRIVATE simplepir)

add_executable(MLWEtoLWE-test test-mlwe2lwe.cpp)
target_compile_options(MLWEtoLWE-test PRIVATE -std=c++11 -O3)
target_link_libraries(MLWEtoLWE-test PRIVATE simplepir)

add_executable(MLWE-test test-mlwe.cpp)
target_compile_options(MLWE-test PRIVATE -std=c++11 -O3)
target_link_libraries(MLWE-test PRIVATE simplepir)

add_executable(mult-test test-mult.cpp)
target_compile_options(mult-test PRIVATE -std=c++11 -O3)
target_link_libraries(mult-test PRIVATE simplepir)
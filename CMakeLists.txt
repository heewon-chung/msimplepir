cmake_minimum_required(VERSION 3.1)

#######################
# project description #
#######################
project(simplepir LANGUAGES CXX)

#######################
# set CMAKE_CXX_FLAGS #
#######################

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
#add_definitions(-D__DEBUG)

####################
# Add Cmake Module #
####################

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

######################
# Build library #
######################

add_subdirectory(src)

#####################
# Make EXE for test #
#####################

add_subdirectory(tests)
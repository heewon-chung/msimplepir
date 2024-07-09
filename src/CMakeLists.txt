
#######################
# headers for library #
#######################
                
# set(HEADER_LIST 
#     "../include/database.h"
#     "../include/util.h"
#     "../include/simplepir.h"
#     "../include/vlhe.h"
# )

#######################
# sources for library #
#######################

set(SRC_LIST 
    "database.cpp"
    "parameter.cpp"
    "util.cpp"
    "simplepir.cpp"
    )

#########################
# make a static library #
#########################

# add_library(simplepir STATIC ${SRC_LIST} ${HEADER_LIST})
add_library(simplepir STATIC ${SRC_LIST})

#########################
# ensure this directory #
#########################

target_include_directories(simplepir PUBLIC ../include)

##################
# add dependency #
##################

# target_link_libraries(izkp PRIVATE ${OPENSSL_LIBRARY})

###############
# using C++11 #
###############

target_compile_options(simplepir PUBLIC -std=c++11 -O3)
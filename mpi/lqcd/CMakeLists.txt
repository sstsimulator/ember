Cmake_minimum_required(VERSION 3.17)
project(LQCD_SKELETON)

# Set the MPI direcotry
set(CMAKE_PREFIX_PATH /Users/knteran/openmpi-4.1.1)

find_package(MPI)

if (MPI_FOUND)
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
else (MPI_FOUND)
    message(SEND_ERROR "This application cannot compile without MPI")
endif (MPI_FOUND)
add_executable(LQCD lqcd.c)
target_link_libraries(LQCD PUBLIC ${MPI_LIBRARIES})

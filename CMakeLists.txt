cmake_minimum_required(VERSION 2.8)
project(LTECell-Search)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11") #设置编译选项
# set(CMAKE_C_COMPILER "gcc")
# set(CMAKE_CXX_COMPILER "g++")
# set(CMAKE_CXX_STANDARD 11)

# #debug 模式
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

#select the release build type by default to get optimization flags
# if(NOT CMAKE_BUILD_TYPE)
#    set(CMAKE_BUILD_TYPE "Release")
#    message(STATUS "Build type not specified: defaulting to release.")
# endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "")

set(BOOST_REQUIRED_COMPONENTS
        program_options
        filesystem
        system
        thread
)

find_package(Boost REQUIRED ${BOOST_REQUIRED_COMPONENTS})

find_package(UHD "3.10.0" REQUIRED)

include_directories(
        ${Boost_INCLUDE_DIRS}
        ${UHD_INCLUDE_DIRS}
)

link_directories(${Boost_LIBRARY_DIRS})

set(SOURCE_FILES main.cpp csync.cpp)

add_executable(LTECell-Search ${SOURCE_FILES})

# set(CTOCMAKE_RUNTIME_OUTPUT_DIRERY ${CMAKE_SOURCE_DIR})

target_link_libraries(LTECell-Search ${UHD_LIBRARIES} ${Boost_LIBRARIES} pthread)

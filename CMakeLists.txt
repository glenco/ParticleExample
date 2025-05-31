cmake_minimum_required(VERSION 3.0...3.8)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

####
# project settings
####

project(TestParicles)

set(sample_VERSION_MAJOR 1)
set(sample_VERSION_MINOR 0)


if(CMAKE_BUILD_TYPE MATCHES Debug)
	message(STATUS "Debug build")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif()
if(CMAKE_BUILD_TYPE MATCHES Release)
	message(STATUS "Release build")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Os -DNDEBUG")
endif()

####
# headers
####

set(HEADERS
)


####
# sources
####

set(SOURCES
	testparticles.cpp
)


####
# dependencies
####

find_package(GLAMER NO_MODULE REQUIRED)

include_directories(${GLAMER_INCLUDE_DIRS})


####
# build
####

add_executable(TestParticles ${HEADERS} ${SOURCES})

target_link_libraries(TestParticles ${GLAMER_LIBRARIES})



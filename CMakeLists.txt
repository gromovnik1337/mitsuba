CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(Mitsuba)

file(GLOB_RECURSE SRCS *.cpp)
file(GLOB_RECURSE HDRS *.h)

include_directories(include)
include_directories(src)

ADD_EXECUTABLE(all ${SRCS} ${HDRS})
TARGET_LINK_LIBRARIES(all)

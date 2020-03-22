cmake_minimum_required(VERSION 3.10)
project(demo)
set(CMAKE_CXX_STANDARD 11)
find_package(OpenCV REQUIRED)

include_directories(include)

file(GLOB SOURCES "src/*.cpp")

add_executable(myexec src/levelset.cpp src/main.cpp)


target_link_libraries( myexec ${OpenCV_LIBS} )
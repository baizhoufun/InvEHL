cmake_minimum_required (VERSION 3.0)
project(test)

option(iniReader_test "abc" OFF)
option(eikonal_test "abc" OFF)
option(pde_test "abc" ON)

add_compile_options(-std=c++11)
add_compile_options(-O3)

if(iniReader_test)
message("iniReader_test is ${iniReader_test}")
add_library(iniReaderLib1 io/iniReader.cpp)
add_library(iniReaderLib2 io/ini.c)
add_executable(main test/iniReader_test.cpp)
target_link_libraries(main PRIVATE iniReaderLib1 iniReaderLib2)
else()
message("iniReader_test is ${iniReader_test}")
endif()

if(eikonal_test)
message("eikonal_test is ${eikonal_test}")
add_library(iniReaderLib1 io/iniReader.cpp)
add_library(iniReaderLib2 io/ini.c)
find_package(OpenCV REQUIRED)
add_library(eikonalLib eikonal/eikonal.cpp)
add_executable(main test/eikonal_test.cpp)
target_link_libraries(main PRIVATE iniReaderLib1 iniReaderLib2 eikonalLib ${OpenCV_LIBS})
else()
message("eikonal_test is ${eikonal_test}")
endif()

if(pde_test)

message("pde_test is ${pde_test}")

find_package (Threads)



find_package(Eigen3 3.1.0 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

add_library(ioLib io/iniReader.cpp  io/ini.c io/ioEigen.cpp)

find_package(OpenCV REQUIRED)
add_library(eikonalLib eikonal/eikonal.cpp)


add_library(pdeLib1 pde/mesh.cpp)

add_executable(main test/pde_test.cpp)
target_link_libraries(main PRIVATE ioLib eikonalLib pdeLib1 ${OpenCV_LIBS} ${EIGEN3_LIBS} ${CMAKE_THREAD_LIBS_INIT})

else()

message("pde_test is ${pde_test}")

endif()
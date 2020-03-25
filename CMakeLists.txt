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

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
message("OPENMP FOUND")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

endif()



find_package(Eigen3 3.1.0 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

set(ROOT ${PROJECT_SOURCE_DIR})
set(SRC "${ROOT}/src")
include_directories("${SRC}")

add_library(ioLib 
${SRC}/io/iniReader.cpp  
${SRC}/io/ini.c 
${SRC}/io/ioEigen.cpp 
${SRC}/io/utilities.cpp)

find_package(OpenCV REQUIRED)
add_library(eikonalLib ${SRC}/eikonal/eikonal.cpp)

add_library(pdeLib
${SRC}/io/ioEigen.cpp 
${SRC}/pde/function.cpp
${SRC}/pde/mesh.cpp 
${SRC}/pde/data.cpp
${SRC}/pde/tfe.cpp)

add_executable(main ${ROOT}/test/pde_test.cpp)
target_link_libraries(main PRIVATE ioLib eikonalLib pdeLib ${OpenCV_LIBS} ${EIGEN3_LIBS} ${CMAKE_THREAD_LIBS_INIT})

else()

endif()
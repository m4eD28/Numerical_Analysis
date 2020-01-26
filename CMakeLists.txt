cmake_minimum_required(VERSION 3.13)

project(Numerical_Analysis CXX)
# find_package(PythonLibs 2.7)
# include_directories(${PYTHON_INCLUDE_DIRS})
set(CMAKE_CXX_STANDARD 11)
# set(CMAKE_CXX_FLAGS "-o3 -std=c++11 -lpython2.7")
set(CMAKE_CXX_FLAGS "-o3 -std=c++11")
# set(PYTHONINTERP_FOUND)
# set(PYTHON_EXECUTABLE)

add_subdirectory(src)
add_subdirectory(exe)

# add_executable(a.out -lstdc++ -lpython2.7 -std=c++11 main2.cpp src/LinearAlgebra.cpp src/Algo.cpp )
add_executable(Numerical_Analysis exe/report_winter_1.cpp src/LinearAlgebra.cpp src/Algo.cpp)

# target_compile_options(Numerical_Analysis PUBLIC -lstdc++ -lpython2.7)
# target_link_libraries(Numerical_Analysis PRIVATE ${PYTHON_LIBRARIES})
# target_include_directories(Numerical_Analysis PUBLIC src)
# target_include_directories(Numerical_Analysis PUBLIC ${PYTHON_INCLUDE_DIRS})

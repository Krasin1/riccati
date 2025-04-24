cmake_minimum_required(VERSION 3.10)
project(RiccatiSolver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(solver main.cpp utils.cpp solver.cpp)

# Добавляем флаг оптимизации 
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    message(STATUS "OpenMP not found. Parallelization disabled.")
endif()

# Ищем библиотеку Eigen3, если Eigen не найден, загружаем его с помощью FetchContent
find_package(Eigen3 QUIET)
if(NOT Eigen3_FOUND)
    message(STATUS "Eigen not found in system. Downloading Eigen...")

    include(FetchContent)
    FetchContent_Declare(
        Eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0  
    )
    FetchContent_MakeAvailable(Eigen)
else()
    message(STATUS "Eigen found in system. Using installed version.")
endif()

#CLI11
find_package(CLI11 QUIET)
if(NOT CLI11_FOUND)
    message(STATUS "CLI11 not found in system. Downloading CLI11...")
    include(FetchContent)
    FetchContent_Declare(
        CLI11
        GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
        GIT_TAG v2.5.0
    )
    FetchContent_MakeAvailable(CLI11)
    target_link_libraries(solver PUBLIC CLI11::CLI11)
else()
    message(STATUS "CLI11 found in system. Using installed version.")
endif()

message(STATUS "Eigen include dir: ${EIGEN3_INCLUDE_DIR}")


target_link_libraries(solver PUBLIC Eigen3::Eigen CLI11::CLI11)

cmake_minimum_required(VERSION 3.10)
project(RiccatiSolver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

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

message(STATUS "Eigen include dir: ${EIGEN3_INCLUDE_DIR}")

add_executable(solver
    main.cpp
    functions/func.cpp
    methods/solver.cpp
    methods/runge/runge.cpp
    methods/adams/adams.cpp
    methods/felberg/felberg.cpp
    methods/hemming/hemming.cpp
    methods/inglend/inglend.cpp
    methods/milna/milna.cpp
    methods/nystrom/nystrom.cpp
)

target_link_libraries(solver PUBLIC Eigen3::Eigen)

cmake_minimum_required(VERSION 3.16)

project(NRRRC)

project(
    ${PROJECT_NAME}
    VERSION 2.0.0
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(matplotplusplus
    GIT_REPOSITORY https://github.com/alandefreitas/matplotplusplus
    GIT_TAG origin/master)
FetchContent_MakeAvailable(matplotplusplus)

FetchContent_GetProperties(matplotplusplus)
if(NOT matplotplusplus_POPULATED)
    FetchContent_Populate(matplotplusplus)
    add_subdirectory(${matplotplusplus_SOURCE_DIR} ${matplotplusplus_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

FetchContent_Declare(eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen
    GIT_TAG 3.4.0)
FetchContent_MakeAvailable(eigen)

include_directories(${PROJECT_SOURCE_DIR}/build/_deps/eigen-src/)

add_subdirectory(src)


# target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_20)

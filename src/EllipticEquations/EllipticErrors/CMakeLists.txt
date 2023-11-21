project(EllipticErrors)

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(
    ${PROJECT_NAME}
    INTERFACE
    Eigen3::Eigen
)
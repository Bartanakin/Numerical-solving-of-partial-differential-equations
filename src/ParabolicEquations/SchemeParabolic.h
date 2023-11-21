#pragma once
#include "SetupParabolic.h"
#include <Eigen/Dense>
#include <memory>
#include <ranges>

using TypedSetupParabolic = SetupParabolic<FType>;
using Vector = typename Eigen::VectorX<FType>;
using Matrix = typename Eigen::MatrixX<FType>;
using SetupMatrices = typename std::pair<Matrix, Vector>;

SetupMatrices create_matrices(const TypedSetupParabolic& setup) {
    auto u = Vector(setup.x_num);
    u[0] = setup.boundary_condition;
    u[u.size() - 1] = setup.boundary_condition;
    for (int i = 1; i < u.size() - 1; i++) {
        u[i] = setup.f(setup.dx * i, 0);
    }

    Matrix A = std::move(Matrix::Identity(setup.x_num - 2, setup.x_num - 2));
    auto dx_square = setup.dx * setup.dx;
    auto p1 = -(2 * setup.alpha / dx_square + 1 / setup.dt);
    auto p2 = setup.alpha / dx_square;
    for (int i = 0; i < setup.x_num - 2; i++) {
        A(i, i) = p1;
        if (i != 0) {
            A(i - 1, i) = p2;
        } 

        if (i != setup.x_num - 3) {
            A(i + 1, i) = p2;
        } 
    }

    return {std::move(A), std::move(u)};
}

Vector create_RH(const Vector& u, const TypedSetupParabolic& setup) {
    auto d = Vector(setup.x_num - 2);

    auto dx_square = setup.dx * setup.dx;
    auto p1 = 2 * (1- setup.alpha) / dx_square - 1 / setup.dt;
    auto p2 = (setup.alpha - 1) / dx_square;

    for (int i = 1; i < u.size() - 1; i++) {
        d[i - 1] = p1 * u[i] + p2 * (u[i - 1] + u[i + 1]);
    }

    return d;
}

std::pair<Matrix, Matrix> scheme_parabolic(const TypedSetupParabolic& setup) {
    auto u_matrix = Matrix(setup.t_out_num, setup.x_num);
    auto u_matrix_exact = Matrix(setup.t_out_num, setup.x_num);
    auto [A, u] = std::move(create_matrices(setup));
    // u_matrix.row(0) = u;
    for (int i = 0; i < setup.x_num; i++) {
        u_matrix_exact(0, i) = setup.exact(i * setup.dx, 0);
    }

    for (int t = 1; t < setup.t_num; t++) {
        auto d = create_RH(u, setup);
        Vector next_d = A.colPivHouseholderQr().solve(d);
        for (int i = 1; i < u.size() - 1; i++) {
            u[i] = next_d[i - 1];
        }

        if ((t - 1) % setup.output_freq == 0) {
            for (int i = 0; i < setup.x_num; i++) {
                u_matrix_exact((t - 1) / setup.output_freq, i) = setup.exact(i*setup.dx, t*setup.dt);
            }
            // std::cout << t << std::endl;
            u_matrix.row((t - 1) / setup.output_freq) = u;
        }
    }

    return {std::move(u_matrix), std::move(u_matrix_exact)};
}
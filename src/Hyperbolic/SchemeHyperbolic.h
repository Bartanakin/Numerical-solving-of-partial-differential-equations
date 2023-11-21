#pragma once
#include "SetupHyperbolic.h"
#include <Eigen/Dense>
#include <memory>
#include <ranges>

using TypedSetupHyperbolic = SetupHyperbolic<FType>;
using Vector = typename Eigen::VectorX<FType>;
using Matrix = typename Eigen::MatrixX<FType>;
using SetupMatrices = typename std::pair<Matrix, Vector>;

SetupMatrices create_matrices(const TypedSetupHyperbolic& setup) {
    auto u = Vector(setup.x_num);
    Matrix A = Matrix();

    return {std::move(A), std::move(u)};
}

Vector create_RH(const Vector& u, const Vector& u_prev, const TypedSetupHyperbolic& setup) {
    auto d = Vector(setup.x_num - 2);
    auto r = setup.dt * setup.dt / (setup.dx * setup.dx);

    return d;
}

FType fictional(FType x, const TypedSetupHyperbolic& setup) {
    return x;
}

std::pair<Matrix, Matrix> scheme_parabolic(const TypedSetupHyperbolic& setup) {
    auto u_matrix = Matrix(setup.t_out_num, setup.x_num);
    auto u_matrix_exact = Matrix(setup.t_out_num, setup.x_num);
    auto [A, u] = std::move(create_matrices(setup));
    Vector u_prev = std::move(Vector::Zero(u.size()));
    for (int i = 0; i < setup.x_num; i++) {
        auto xi = setup.t_range.first + setup.dx * i;

        u_prev[i] = fictional(xi, setup);
    }

    for (int t = 1; t < setup.t_num; t++) {
        auto d = create_RH(u, u_prev, setup);

        if ((t - 1) % setup.output_freq == 0) {
            u_matrix.row((t - 1) / setup.output_freq) = u;
        }
    }

    return {std::move(u_matrix), std::move(u_matrix_exact)};
}
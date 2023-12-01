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
    u[0] = setup.boundary_condition.first(0);
    u[u.size() - 1] = setup.boundary_condition.second(0);
    for (int i = 1; i < u.size() - 1; i++) {
        u[i] = setup.initial(setup.dx * i);
    }

    Matrix A = std::move(Matrix::Identity(setup.x_num - 2, setup.x_num -2));
    auto r = setup.dt * setup.dt / setup.dx / setup.dx;
    for (int i = 1; i < u.size() - 1; i++) {
        if (i != 1) {
            A(i - 2, i - 1) = - r * setup.alpha;
        }

        A(i - 1, i - 1) = r * 2.f * setup.alpha + 1.f;
        
        if (i != u.size() - 2) {
            A(i, i - 1) = - r * setup.alpha;
        }
    }

    return {std::move(A), std::move(u)};
}

Vector create_RH(const Vector& u, const Vector& u_prev, const TypedSetupHyperbolic& setup) {
    auto d = Vector(setup.x_num - 2);
    auto r = setup.dt * setup.dt / setup.dx / setup.dx;

    // for (int i = 1; i < u.size() - 1; i++) {
    //     d[i - 1] = r * (u[i + 1] - 2* u[i] + u[i - 1]) + 2 * u[i] - u_prev[i];
    // }
    for (int i = 1; i < u.size() - 1; i++) {
        d[i - 1] = (1 - 2.f * setup.alpha) * r * (u[i + 1] - 2* u[i] + u[i - 1])
                    + setup.alpha * r * (u_prev[i + 1] - 2* u_prev[i] + u_prev[i - 1])
                    + 2 * u[i] - u_prev[i];
    }



    return d;
}

FType fictional(FType x, const TypedSetupHyperbolic& setup) {
    auto Lhp = (1 - 2.f * setup.alpha) * (setup.initial(x - setup.dx) - 2 * setup.initial(x) + setup.initial(x + setup.dx));
    Lhp += setup.alpha * (setup.initial_der(x - setup.dx) - 2 * setup.initial_der(x) + setup.initial_der(x + setup.dx));
    Lhp *= 1 / (setup.dx * setup.dx);

    return setup.initial(x) - setup.dt * setup.initial_der(x) + setup.dt * setup.dt * 0.5f *Lhp;
}

std::pair<Matrix, Matrix> scheme_parabolic(const TypedSetupHyperbolic& setup) {
    auto u_matrix = Matrix(setup.t_out_num, setup.x_num);
    auto u_matrix_exact = Matrix(setup.t_out_num, setup.x_num);
    auto [A, u] = std::move(create_matrices(setup));
    for (int i = 0; i < setup.x_num; i++) {
        u_matrix_exact(0, i) = setup.exact(i * setup.dx, 0);
    }

    Vector u_prev = std::move(Vector::Zero(u.size()));
    for (int i = 0; i < setup.x_num; i++) {
        auto xi = setup.t_range.first + setup.dx * i;

        u_prev[i] = fictional(xi, setup);
    }

    for (int t = 1; t < setup.t_num; t++) {
        auto d = create_RH(u, u_prev, setup);
        Vector next_d = A.colPivHouseholderQr().solve(d);

        for (int i = 0; i < u.size(); i++) {
            u_prev[i] = u[i];
        }

        u[0] = setup.boundary_condition.first(t * setup.dt);
        u[u.size() - 1] = setup.boundary_condition.second(t * setup.dt);
        for (int i = 1; i < u.size() - 1; i++) {
            u[i] = next_d[i - 1];
        }

        if ((t - 1) % setup.output_freq == 0) {
            for (int i = 0; i < setup.x_num; i++) {
                u_matrix_exact((t - 1) / setup.output_freq, i) = setup.exact(i*setup.dx, t*setup.dt);
            }

            u_matrix.row((t - 1) / setup.output_freq) = u;
        }
    }

    return {std::move(u_matrix), std::move(u_matrix_exact)};
}
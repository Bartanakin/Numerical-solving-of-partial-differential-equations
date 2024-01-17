#pragma once
#include "SetupElliptic.h"
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>

template<size_t N>
bool isIndexForInitial(size_t i, const SetupElliptic<FType>& setup) {
    if (i == 0 || i == setup.x_num - 1 || i == N - setup.x_num  || i == N - 1) {
        return true;
    } else if (i < setup.x_num || i >= N - setup.x_num) {
        return true;
    } else if (i % setup.y_num == 0 || i % setup.y_num == setup.y_num - 1) {
        return true;
    }

    return false;
}

template<size_t N, typename DoubleIterator> 
Eigen::VectorX<FType> scheme_elliptic(const SetupElliptic<FType>& setup) {
    using VecN = Eigen::VectorX<FType>;
    using MatN = Eigen::SparseMatrix<FType>;
    const auto size = N;

    auto hx = - 1. / (setup.dx * setup.dx); 
    auto hy = - 1. / (setup.dy * setup.dy);
    auto r = - 2. * (hx + hy);

    auto b = VecN(N);

    for (DoubleIterator it; it.isValid(); ++it) {
        auto boundary = setup.boundaryCondition(setup.dx * it.i1(), setup.dy * it.i2());
        if (isIndexForInitial<N>(it.total(), setup)) { // Nodes on the bound
            b(it.total()) = boundary;
        } else {
            b[it.total()] = setup.f(it.i1() * setup.dx, it.i2() * setup.dy);// Nodes inside
            if (it.total() < setup.x_num * 2 || it.total() >= N - 2 * setup.x_num) { // Nodes right next to horizontal bounds
                b(it.total()) -= boundary * hx;
            }

            if (it.total() % setup.y_num == 1 || it.total() % setup.y_num == setup.y_num - 2) { // Nodes right next to vertical bounds
                b(it.total()) -= boundary * hy;
            }
        }

    }

    std::vector<Eigen::Triplet<FType>> tripletList;
    tripletList.reserve(5 * N - 2 * 4 * setup.x_num - 2 * 4 * setup.y_num + 4 * 4);
    for (int i = 0; i < N; i++) {
        if (isIndexForInitial<N>(i, setup)) { // Nodes on the bound
            tripletList.emplace_back(i, i, 1.);
        } else {
            tripletList.emplace_back(i, i, r);
            if (i < N - 1 && !isIndexForInitial<N>(i + 1, setup)) { // Two diagonals next to the main diagonal
                tripletList.emplace_back(i + 1, i, hx);
                tripletList.emplace_back(i, i + 1, hx);
            }

            if (i < N - setup.x_num && !isIndexForInitial<N>(i + setup.x_num, setup)) { // partial diagonals
                tripletList.emplace_back(i + setup.x_num, i, hy);
            }

            if (i < N - setup.y_num && !isIndexForInitial<N>(i + setup.y_num, setup)) { // partial diagonals
                tripletList.emplace_back(i, i + setup.y_num, hy);
            }
        }
    }

    auto A = MatN(N, N);
    A.setFromTriplets(tripletList.begin(), tripletList.end());

    Eigen::SimplicialCholesky<MatN> chol(A);
    Eigen::VectorX<FType> x = chol.solve(b);
    // std::cout << A;
    // std::cout << b;

    return x;
}
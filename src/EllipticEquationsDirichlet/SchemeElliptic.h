#pragma once
#include "SetupElliptic.h"
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>

template<size_t N>
bool isIndexForInitial(size_t i, const SetupElliptic<FType>& setup) {
    if (i == N - 1) {
        return false;
    } else if (i == 0 || i == setup.x_num - 1 || i == N - setup.x_num) {
        return true;
    } else if (i < setup.x_num) {
        return true;
    } else if (i % setup.y_num == 0) {
        return true;
    }

    return false;
}
template<size_t N>
bool isIndexForDirichlet(size_t i, const SetupElliptic<FType>& setup) {
    if (i == N - 1) {
        return true;
    } else if (i == 0 || i == setup.x_num - 1 || i == N - setup.x_num) {
        return false;
    } else if (i >= N - setup.x_num) {
        return true;
    } else if (i % setup.y_num == setup.y_num - 1) {
        return true;
    }

    return false;
}

template<size_t N>
bool isIndexForXDirichlet(size_t i, const SetupElliptic<FType>& setup) {
    if (i == N - 1) {
        return true;
    } else if (i == 0 || i == setup.x_num - 1 || i == N - setup.x_num) {
        return false;
    } else if (i >= N - setup.x_num) {
        return false;
    } else if (i % setup.y_num == setup.y_num - 1) {
        return true;
    }

    return false;
}

template<size_t N>
bool isIndexForYDirichlet(size_t i, const SetupElliptic<FType>& setup) {
    if (i == N - 1) {
        return true;
    } else if (i == 0 || i == setup.x_num - 1 || i == N - setup.x_num) {
        return false;
    } else if (i >= N - setup.x_num) {
        return true;
    } else if (i % setup.y_num == setup.y_num - 1) {
        return false;
    }

    return false;
}

template<size_t N>
bool isIndexForBothDirichlet(size_t i, const SetupElliptic<FType>& setup) {
    if (i == N - 1) {
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
        auto f = setup.f(it.i1() * setup.dx, it.i2() * setup.dy);
        auto boundary = setup.boundaryCondition(setup.dx * it.i1(), setup.dy * it.i2());
        auto dirichlet = 1;
        if (isIndexForInitial<N>(it.total(), setup)) { // Nodes on the bound
            b(it.total()) = boundary;

            continue;
        }
        
        b[it.total()] = f;// Nodes inside
        if (isIndexForBothDirichlet<N>(it.total(), setup)) { // Nodes on the bound
            b(it.total()) = 2 * (1. / setup.dx + 1. / setup.dy) * dirichlet;
        } else if (isIndexForXDirichlet<N>(it.total(), setup)) {
            b(it.total()) = 2 * dirichlet / setup.dx;
        } else if (isIndexForYDirichlet<N>(it.total(), setup)) {
            b(it.total()) = 2 * dirichlet / setup.dy;
        } else {
            if (it.total() < setup.x_num * 2) { // Nodes right next to the lower horizontal bound
                b(it.total()) -= boundary * hy;
            }

            if (it.total() % setup.y_num == 1) { // Nodes right next to the left vertical bound
                b(it.total()) -= boundary * hx;
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
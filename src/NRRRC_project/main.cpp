#include <iostream>
#include <functional>
#include <vector>
#include <Eigen/Dense>
#include <filesystem>
#include <fstream>
#include <iterator>

using Solution = std::function<Eigen::VectorXd(double, const Eigen::VectorXd&)>;

constexpr const size_t steps = 50 * 40 + 1;
//constexpr const size_t steps = 100 + 1;
constexpr const double step_size = 0.0005;
constexpr const size_t NODE_COUNT = 41;
constexpr const size_t ELEMENT_COUNT = NODE_COUNT - 2;
constexpr const double k = 1.;
constexpr const double h = 1. / (static_cast<double>(NODE_COUNT - 1));
constexpr const size_t T_OUTPUT_INTERVAL = 40;
constexpr const double X_OUTPUT_STEP = 0.025;
constexpr const size_t T_OUTPUT_SIZE = (steps / T_OUTPUT_INTERVAL) + 1;
constexpr const size_t X_OUTPUT_SIZE = static_cast<size_t>(1. / X_OUTPUT_STEP) + 1;
constexpr const double T_OUTPUT_STEP = (step_size * T_OUTPUT_INTERVAL);

Eigen::MatrixXd rk4(
    const Solution& f,
    const Eigen::VectorXd& t,
    const Eigen::VectorXd& x0
) {
    auto x = Eigen::MatrixXd(t.size(), x0.size());
    x.row(0).transpose() = x0;
    for (int i = 1; i < t.size(); i++) {
        auto dt = t[i] - t[i - 1];
        Eigen::VectorXd xi = x.row(i - 1).transpose();
        Eigen::VectorXd k1 = dt * f(t[i - 1], xi);
        assert(k1.size() == xi.size());
        Eigen::VectorXd k2 = dt * f(t[i - 1] + 0.5 * dt, xi + 0.5 * k1);
        Eigen::VectorXd k3 = dt * f(t[i - 1] + 0.5 * dt, xi + 0.5 * k2);
        Eigen::VectorXd k4 = dt * f(t[i - 1] * dt, xi + k3);

        auto sol = xi + (k1 + 2. * k2 + 2. * k3 + k4) / 6.;
        x.row(i) = sol.transpose();
    }

    return x.transpose();
}

double exact(double t, double x) {
    return std::exp(-M_PI * M_PI * t) * std::sin(M_PI * x);
}

Eigen::MatrixXd open_euler(
        const Solution& f,
        const Eigen::VectorXd& t,
        const Eigen::VectorXd& x0
) {
    auto x = Eigen::MatrixXd(t.size(), x0.size());
    x.row(0).transpose() = x0;
    for (int i = 1; i < t.size(); i++) {
        auto dt = t[i] - t[i - 1];
        Eigen::VectorXd xi = x.row(i - 1).transpose();

        x.row(i) = (xi + dt * f(t[i - 1], xi)).transpose();
    }

    return x.transpose();
}

Eigen::MatrixXd closed_euler(
        const Solution& f,
        const Eigen::VectorXd& t,
        const Eigen::VectorXd& x0
) {
    auto x = Eigen::MatrixXd(t.size(), x0.size());
    x.row(0).transpose() = x0;
    for (int i = 1; i < t.size(); i++) {
        auto dt = t[i] - t[i - 1];
        Eigen::VectorXd xi = x.row(i - 1).transpose();

        x.row(i) = f(dt, xi).transpose();
    }

    return x.transpose();
}

double base_func(double x, int i) {
    i++;
    if (x < h * (static_cast<double>(i - 1)) || x > h * (static_cast<double>(i + 1))) {
        return 0.;
    }

    if (x < h * static_cast<double>(i)) {
        return x / h + static_cast<double>(1 - i);
    }

    return -x / h + static_cast<double>(1 + i);
}

double u(const Eigen::VectorXd& phi, double x) {
    assert(phi.size() == ELEMENT_COUNT);
    double sum = 0.;
    // TODO optimise
    for (int j = 0; j < phi.size(); j++) {
        sum += phi[j] * base_func(x, j);
    }

    return sum;
}

int main() {

    auto t = Eigen::VectorXd(steps);
    t[0] = 0.f;
    for (int i = 1; i < steps; i++) {
        t[i] = t[i - 1] + step_size;
    }

    auto A = Eigen::MatrixXd(ELEMENT_COUNT, ELEMENT_COUNT);
    A.fill(0.);
    for (int i = 0; i < ELEMENT_COUNT; i++) {
        A(i, i) = 4.;
        if (i != 0) {
            A(i - 1, i) = 1.;
        }

        if (i != ELEMENT_COUNT - 1) {
            A(i + 1, i) = 1.;
        }
    }

    auto B = Eigen::MatrixXd(ELEMENT_COUNT, ELEMENT_COUNT);
    B.fill(0.);
    for (int i = 0; i < ELEMENT_COUNT; i++) {
        B(i, i) = 2.;
        if (i != 0) {
            B(i - 1, i) = -1.;
        }

        if (i != ELEMENT_COUNT - 1) {
            B(i + 1, i) = -1.;
        }
    }

    auto A_inv = A.inverse();
    auto r = 6. / (h * h);
    auto D_QR = A.colPivHouseholderQr();
    Solution solution_func = [&](double t, const Eigen::VectorXd& x) {
        Eigen::VectorXd xd = -r * D_QR.solve(B * x);

        return xd;
    };

    Eigen::MatrixXd J = (6.* step_size / (h * h)) * A_inv * B + Eigen::MatrixXd::Identity(A.rows(), A.cols());
    auto J_QR = J.colPivHouseholderQr();
    Solution closed_solution_func = [&](double t, const Eigen::VectorXd& x) {
        Eigen::VectorXd xd = J_QR.solve(x);

        return xd;
    };

    auto x0 = Eigen::VectorXd(ELEMENT_COUNT);
    for (int i = 1; i < NODE_COUNT - 1; i++) {
        auto x = i * h;
        x0[i - 1] = std::sin(M_PI * x);
    }

//    auto solution = open_euler(solution_func, t, x0);
//    auto solution = rk4(solution_func, t, x0);
    auto solution = closed_euler(closed_solution_func, t, x0);
    auto values = Eigen::MatrixXd(NODE_COUNT - 1, solution.cols());

    // output
    std::filesystem::path output_path = "ChartData";
    output_path /= "project_NRRRC";
    std::array<std::ofstream, 3> files = {
            std::ofstream((output_path / "numeric_result.txt").string()),
            std::ofstream((output_path / "exact.txt").string()),
            std::ofstream((output_path / "error.txt").string())
    };
    for (auto& file : files) {
        file << T_OUTPUT_SIZE << "\n";
        file << X_OUTPUT_SIZE << "\n";
    }
    std::vector<double> axis(T_OUTPUT_SIZE);
    std::vector<double> axis2(X_OUTPUT_SIZE);
    for (auto& file : files) {
        std::generate(axis.begin(), axis.end(), [beg = -T_OUTPUT_STEP, step = T_OUTPUT_STEP] () mutable { return beg += step; });
        std::generate(axis2.begin(), axis2.end(), [beg = -X_OUTPUT_STEP, step = X_OUTPUT_STEP] () mutable { return beg += step; });
        std::copy(axis.begin(), axis.end(), std::ostream_iterator<double>(file, " "));
        std::copy(axis2.begin(), axis2.end(), std::ostream_iterator<double>(file, " "));
    }

    auto maxError = 0.;
    for (int i = 0; i < T_OUTPUT_SIZE; i++) {
        for (int j = 0; j < X_OUTPUT_SIZE; j++) {
            auto numeric = u(solution.col(i * T_OUTPUT_INTERVAL), j * X_OUTPUT_STEP);
            auto exact_s = exact(i * step_size * T_OUTPUT_INTERVAL, j * X_OUTPUT_STEP);

            files[0] << numeric << " ";
            files[1] << exact_s << " ";
            files[2] << std::abs(numeric - exact_s) << " ";

            if (std::abs(numeric - exact_s) > maxError) {
                maxError = std::abs(numeric - exact_s);
            }
        }
    }

    std::cout << "Max error = " << maxError << std::endl;
}

#include <iostream>
#include <fstream>
#include <filesystem>
#include <algorithm>
#include <iterator>
#include <chrono>
#include "../../Presets.h"
#include "../SetupElliptic.h"
#include "../../DoubleLoopIterator.h"
#include "../SchemeElliptic.h"

int main() {
    constexpr const SetupElliptic<FType> setup = {0.01f};
    // std::cout << SetupElliptic<FType>::exact(.5f, .5f);

    auto time_point = std::chrono::steady_clock::now();
    constexpr const auto size = setup.x_num * setup.y_num;
    auto x = scheme_elliptic<size, DoubleLoopIterator<setup.x_num, setup.y_num>>(setup);
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - time_point) << "\n";
    Eigen::MatrixX<FType> grid(setup.x_num, setup.y_num);
    for (DoubleLoopIterator<setup.x_num, setup.y_num> it; it.isValid(); ++it) {
        grid(it.i1(), it.i2()) = x(it.total());
    }

    // output
    std::filesystem::path output_path = "ChartData";
    output_path /= "Elliptic";
    std::array<std::ofstream, 3> files = {
        std::ofstream((output_path / "numeric_result.txt").string()),
        std::ofstream((output_path / "exact.txt").string()),
        std::ofstream((output_path / "error.txt").string())
    };
    for (auto& file : files) {
        file << setup.x_num << "\n";
        file << setup.y_num << "\n";
    }
    std::vector<FType> axis(setup.x_num);
    for (auto& file : files) {
        std::generate(axis.begin(), axis.end(), [beg = setup.x_range.first - setup.dx, step = setup.dx] () mutable { return beg += step; });
        for (auto x : {0, 1}) {
            std::copy(axis.begin(), axis.end(), std::ostream_iterator<FType>(file, " "));
            file << "\n";
        }
    }
    for (DoubleLoopIterator<setup.x_num, setup.y_num> it; it.isValid(); ++it) {
        files[0] << grid(it.i1(), it.i2()) << " ";
        files[1] << setup.exact(axis[it.i1()], axis[it.i2()]) << " ";
        files[2] << abs(grid(it.i1(), it.i2()) - setup.exact(axis[it.i1()], axis[it.i2()])) << " ";
    }


    return 0;
}
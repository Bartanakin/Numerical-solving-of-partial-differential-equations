#include <iostream>
#include <fstream>
#include <filesystem>
#include <algorithm>
#include <iterator>
#include <chrono>
#include "../../Presets.h"
#include "../../DDoubleLoopIterator.h"
#include "../SchemeHyperbolic.h"

int main() {
    const TypedSetupHyperbolic setup = {20, 1600, 100, 0.5};

    auto time_point = std::chrono::steady_clock::now();
    auto [u_matrix, u_matrix_exact] = std::move(scheme_parabolic(setup));
    // std::cout << u_matrix << std::endl;
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - time_point) << "\n";
    const auto t_size = setup.t_out_num;

    // output
    std::filesystem::path output_path = "ChartData";
    output_path /= "Hyperbolic";
    std::array<std::ofstream, 3> files = {
        std::ofstream((output_path / "numeric_result.txt").string()),
        std::ofstream((output_path / "exact.txt").string()),
        std::ofstream((output_path / "error.txt").string())
    };
    for (auto& file : files) {
        file << setup.x_num << " " << t_size <<  "\n";
    }

    std::vector<FType> x_axis(setup.x_num);
    for (int i = 0; i < x_axis.size(); i++) {
        x_axis[i] = setup.dx * i;
    }

    std::vector<FType> t_axis(t_size);
    for (int i = 0; i < t_axis.size(); i++) {
        t_axis[i] = setup.dt * setup.output_freq * i;
    }
    
    for (auto& file : files) {
        for (auto space : {&x_axis, &t_axis}) {
            std::copy(space->begin(), space->end(), std::ostream_iterator<FType>(file, " "));
            file << "\n";
        }
    }
    for (auto it = DDoubleLoopIterator<>(t_size, setup.x_num); it.isValid(); ++it) {
        files[0] << u_matrix(it.i2(), it.i1()) << " ";
        files[1] << u_matrix_exact(it.i2(), it.i1()) << " ";
        files[2] << abs(u_matrix(it.i2(), it.i1()) - u_matrix_exact(it.i2(), it.i1())) << " ";
    }


    return 0;
}
#include <iostream>
#include "../../Presets.h"
#include "../SetupElliptic.h"
#include "../../DoubleLoopIterator.h"
#include "../SchemeElliptic.h"
#include <chrono>

using FType = float;

template<FType step>
FType maxError() {
    constexpr const auto setup = SetupElliptic<FType>(step);
    auto time_point = std::chrono::steady_clock::now();
    auto x = scheme_elliptic<setup.size(), DoubleLoopIterator<setup.x_num, setup.y_num>>(setup);
    std::cout << "Time for step " << step << ": " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - time_point) << "\n";

    FType er = 0;
    for (DoubleLoopIterator<setup.x_num, setup.y_num> it; it.isValid(); ++it) {
        auto check_er = abs(x(it.total()) - setup.exact(it.i1() * setup.dx, it.i2() * setup.dy));
        if (check_er > er) 
            er = check_er;
    }

    return er;
}

int main() {
    constexpr const size_t TRIALS = 5;
    constexpr const std::array<FType, TRIALS> steps = { 0.1f, 0.05f, 0.01f, 0.005f, 0.001f };
    std::array<FType, TRIALS> errors;
    errors[0] = maxError<steps[0]>();
    errors[1] = maxError<steps[1]>();
    errors[2] = maxError<steps[2]>();
    errors[3] = maxError<steps[3]>();
    errors[4] = maxError<steps[4]>();
    for (int i = 0; i < TRIALS; i++) {
        std::cout << "Max Error for step " << steps[i] << ": " << errors[i] << "\n";
    }

    return 0;
}
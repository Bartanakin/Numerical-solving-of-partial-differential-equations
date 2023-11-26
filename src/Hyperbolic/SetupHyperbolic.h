#pragma once
#include <math.h>
#include <utility>
#include <functional>
#include "../Presets.h"

#define _USE_MATH_DEFINES = 1

template<typename T>
class SetupHyperbolic {
    public:

    SetupHyperbolic(size_t x_num = 100, size_t t_num = 1000, int output_freq = 10, T alpha = 0) :
        x_range({0, 1}),
        t_range({0, 1}),
        boundary_condition({
            [] (T t) { return 0; },
            [] (T t) { return 0; }
        }),
        dx(static_cast<FType>(1) / static_cast<FType>(x_num)),
        dt(static_cast<FType>(1) / static_cast<FType>(t_num)),
        alpha(alpha),
        x_num(x_num + 1),
        t_num(t_num + 1),
        output_freq(output_freq),
        t_out_num(static_cast<int>(t_num) / output_freq)
    {}

    size_t size() const {
        return this->x_num * this->t_num;
    }

    static T initial(T x) {
        return 0.125 * std::sin(M_PI * x);
    }

    static T initial_der(T x) {
        return 0;
    }

    static T exact(T x, T t) {
        return 0.125 * std::sin(M_PI * x) * std::cos(M_PI * t);
    }

    const std::pair<T, T> x_range;
    const std::pair<T, T> t_range;
    const std::pair<std::function<T(T)>, std::function<T(T)>> boundary_condition;
    const T dx;
    const T dt;
    const T alpha;
    const size_t x_num;
    const size_t t_num;
    const int output_freq;
    const int t_out_num;
};
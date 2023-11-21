#pragma once
#include <math.h>
#include <utility>
#include "../Presets.h"

#define _USE_MATH_DEFINES = 1

template<typename T>
class SetupParabolic {
    public:

    constexpr SetupParabolic(size_t x_num, size_t t_num, int output_freq = 100, T alpha = 0) :
        x_range({0, 1}),
        t_range({0, 1}),
        boundary_condition(0),
        dx(static_cast<FType>(1) / static_cast<FType>(x_num)),
        dt(static_cast<FType>(1) / static_cast<FType>(t_num)),
        alpha(alpha),
        x_num(x_num + 1),
        t_num(t_num + 1),
        output_freq(output_freq),
        t_out_num(static_cast<int>(t_num) / output_freq)
    {}

    constexpr SetupParabolic(T dx) :
        SetupParabolic(dx, dx*dx)
    {}

    constexpr size_t size() const {
        return this->x_num * this->t_num;
    }

    static T f(T x, T t) {
        return std::sin(M_PI * x);
    }

    static T exact(T x, T t) {
        return std::exp(-M_PI * M_PI * t) * std::sin(M_PI * x);
    }

    const std::pair<T, T> x_range;
    const std::pair<T, T> t_range;
    const T boundary_condition;
    const T dx;
    const T dt;
    const T alpha;
    const size_t x_num;
    const size_t t_num;
    const int output_freq;
    const int t_out_num;
};
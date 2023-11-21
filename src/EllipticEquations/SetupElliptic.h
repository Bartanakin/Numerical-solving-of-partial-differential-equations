#pragma once
#include <math.h>
#include <utility>

#define _USE_MATH_DEFINES = 1

template<typename T>
class SetupElliptic {
    public:

    constexpr SetupElliptic(T dx, T dy) :
        x_range({0, 1}),
        y_range({0, 1}),
        boundary_condition(0),
        dx(dx),
        dy(dy),
        x_num(static_cast<size_t>((this->x_range.second - this->x_range.first) / dx) + 1),
        y_num(static_cast<size_t>((this->y_range.second - this->y_range.first) / dy) + 1)
    {}

    constexpr SetupElliptic(T dx) :
        SetupElliptic(dx, dx)
    {}

    constexpr size_t size() const {
        return this->x_num * this->y_num;
    }

    static T f(T x, T y) {
        return 2 * M_PI * M_PI * std::sin(M_PI * x) * std::sin(M_PI * y);
    }

    static T exact(T x, T y) {
        return std::sin(M_PI * x) * std::sin(M_PI * y);
    }

    const std::pair<T, T> x_range;
    const std::pair<T, T> y_range;
    const T boundary_condition;
    const T dx;
    const T dy;
    const size_t x_num;
    const size_t y_num;
};
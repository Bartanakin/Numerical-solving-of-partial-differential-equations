#pragma once

template<size_t D1, size_t D2, typename T = int, T beg = 0>
class DoubleLoopIterator {
public:
    constexpr static T N = D1 * D2;

    DoubleLoopIterator operator++() {
        this->ij++;
        this-> i = this->ij / D1; 
        this-> j = this->ij % D1; 

        return *this;
    }

    inline bool isValid() const {
        return this->ij < N;
    }

    inline T total() const noexcept {
        return this->ij;
    }

    inline T i1() const noexcept {
        return this->i;
    }

    inline T i2() const noexcept {
        return this->j;
    }
private:
    T ij = beg;
    T i = beg;
    T j = beg;
} ;
#pragma once

template<typename T = int>
class DDoubleLoopIterator {
public:
    DDoubleLoopIterator(T D1, T D2, T beg = 0) :
        D1(D1),
        D2(D2),
        N(D1 * D2),
        ij(beg),
        i(beg),
        j(beg)
    {}

    DDoubleLoopIterator operator++() {
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
    const T D1;
    const T D2;
    const T N;
    T ij;
    T i;
    T j;
} ;
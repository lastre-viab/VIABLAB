#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>
#include <string>
#include <algorithm>

using std::string;
using std::max;
using std::min;

const double PLUS_INF = 1e+15;
const double pi = 3.14159265359;
const double numZERO = 0.0000000001;

template <typename T>
T sign(T value) {
    return (value < 0.0) - (value > 0.0);
}

inline bool doubleEquals(double a, double b, double precision=numZERO) {
    return fabs(b - a) < precision;
}

#endif /* CONSTANTS_H */

#define IGNORE_MATH_EXCEPTIONS
#include "../includes/matrix.h"
#include <iostream>

using namespace GeoMath;

int main() {
    mat3<float> m(1, 1, -1, -2, 0, 1, 0, 2, 1);
    auto inv = Inverse(m);
    std::cout << inv << std::endl;
    std::cout << Dot(inv, m) << std::endl;
    return 0;
}

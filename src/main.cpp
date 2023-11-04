#define IGNORE_MATH_EXCEPTIONS
#include "../includes/matrix.h"
#include <iostream>

// using namespace GeoMath;

int main() {
    GeoMath::mat2<float> mat1(
        std::vector<std::vector<float>>({{1, 2}, {3, 4}}));
    GeoMath::mat2x3<float> mat2(
        std::vector<std::vector<float>>({{1, 2, 3}, {4, 5, 6}}));
    std::cout << GeoMath::dot(mat1, mat2) << std::endl;
    return 0;
}

#include "../includes/matrix.h"
#include <iostream>

using namespace GeoMath;

int main() {
    mat2<float> mat = {1, 2, 3, 4};
    std::cout << determinant(mat);
    return 0;
}

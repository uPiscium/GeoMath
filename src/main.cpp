#include <iostream>
#include "../includes/vector.h"

using namespace GeoMath;

int main()
{
	vec3<float> vx(1, 0, 0);
	vec3<float> vy(0, 1, 0);
	std::cout << cross(vx, vy) << std::endl;
	return 0;
}

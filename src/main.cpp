#include <iostream>
#include "../includes/vector.h"

int main()
{
	GeoMath::vec2<double> vec(1, 1.7320508);
	GeoMath::vec3<double> v3 = GeoMath::vec3<double>(vec);
	std::cout << v3 << std::endl;
	return 0;
}

#include <iostream>
#include <initializer_list>
#include <vector>

namespace GeoMath
{
  template <typename T, size_t Row, size_t Column>
  class Matrix
  {
  private:
    using value_type = T;
    using size_type = size_t;

  public:
    Matrix();
    Matrix(const std::vector<T> &data);
    Matrix(const std::initializer_list<T> &data);
    Matrix(const std::vector<Matrix<T, 1, Column>> &data);
    Matrix(const std::vector<Matrix<T, Row, 1>> &data);
    Matrix(const Matrix &data);
  };
}

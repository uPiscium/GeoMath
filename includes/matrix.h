#pragma once
#include <initializer_list>
#include <iostream>
#include <vector>

#include "exceptions.h"
#include "vector.h"

namespace GeoMath {
template <typename T, size_t Row, size_t Column> class MatrixBase {
protected:
  size_t mSize = Row * Column;
  T *mArray = nullptr;
  VectorBase<T, Column> *mRows = nullptr;

private:
  void Allocate() {
    mArray = new T[mSize];
    for (int i = 0; i < mSize; i++) {
      mArray[i] = static_cast<T>(0);
    }

    mRows = new VectorBase<T, Column>[Row];
    for (int i = 0; i < Row; i++) {
      mRows[i] = VectorBase<T, Column>(&mArray[i * Column]);
    }
  }
  void Free() { delete[] mArray; }
  void Flatten(const std::vector<std::vector<T>> &data) {
    for (int i = 0; i < (Row < data.size() ? Row : data.size()); i++) {
      T *ptr = &mArray[i * Column];
      std::memcpy(ptr, &data[i][0],
                  sizeof(T) *
                      (Column < data[i].size() ? Column : data[i].size()));
    }
  }
  void Flatten(const std::initializer_list<std::initializer_list<T>> &data) {
    for (int i = 0; i < (Row < data.size() ? Row : data.size()); i++) {
      T *ptr = &mArray[i * Column];
      std::memcpy(ptr, &data[i][0],
                  sizeof(T) *
                      (Column < data[i].size() ? Column : data[i].size()));
    }
  }

public:
  MatrixBase() { this->Allocate(); }
  MatrixBase(const std::vector<std::vector<T>> &data) {
    this->Allocate();
    this->Flatten(data);
  }
  MatrixBase(const std::initializer_list<std::initializer_list<T>> &data) {
    this->Allocate();
    this->Flatten(data);
  }
  MatrixBase(const std::vector<T> &data) {
    this->Allocate();
    std::memcpy(mArray, &data[0],
                sizeof(T) * (mSize < data.size() ? mSize : data.size()));
  }
  MatrixBase(const std::initializer_list<T> &data) {
    this->Allocate();
    std::memcpy(mArray, &data[0],
                sizeof(T) * (mSize < data.size() ? mSize : data.size()));
  }
  ~MatrixBase() { this->Free(); }

  VectorBase<T, Row> &operator[](const size_t &idx) {
    if (idx < 0 || Row - 1 < idx) {
      _MATH_THROW(IndexOutOfRange("Index is out of range."));
    }
    return mRows[idx];
  }
};
} // namespace GeoMath

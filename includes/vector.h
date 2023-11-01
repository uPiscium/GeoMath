#include <cmath>
#include <initializer_list>
#include <iostream>
#include <vector>

#include "exceptions.h"

namespace GeoMath {
template <typename T, size_t Comp> class VectorBase {
protected:
  size_t mSize = Comp;
  T *mArray = nullptr;
  bool mSelfAlloc = false;

private:
  void Allocate() {
    mArray = new T[mSize];
    for (int i = 0; i < mSize; i++) {
      mArray[i] = static_cast<T>(0);
    }
    mSelfAlloc = true;
  }
  void Free() {
    if (mSelfAlloc) {
      delete[] mArray;
    }
  }

public:
  VectorBase() { this->Allocate(); }
  VectorBase(const std::vector<T> &data) {
    this->Allocate();
    size_t copysize = sizeof(T) * (mSize <= data.size() ? mSize : data.size());
    std::memcpy(mArray, &data[0], copysize);
  }
  VectorBase(const T *data, const size_t &comps) {
    this->Allocate();
    size_t copysize = sizeof(T) * (mSize < comps ? mSize : comps);
    std::memcpy(mArray, data, copysize);
  }
  explicit VectorBase(T *data) : mArray(data) { ; }
  VectorBase(const VectorBase &data) {
    if (data.mSelfAlloc) {
      this->Allocate();
      std::memcpy(mArray, (const T *)data, sizeof(T) * mSize);
    } else {
      mSelfAlloc = false;
      mArray = data.mArray;
    }
  }
  ~VectorBase() { this->Free(); }

  T &operator[](const size_t &idx) {
    if (idx < 0 || mSize - 1 < idx) {
      _MATH_THROW(IndexOutOfRange("Index is out of range."));
    }
    return mArray[idx];
  }
  const T &operator[](const size_t &idx) const {
    if (idx < 0 || mSize - 1 < idx) {
      _MATH_THROW(IndexOutOfRange("Index is out of range."));
    }
    return mArray[idx];
  }

  VectorBase<T, Comp> &operator=(const VectorBase<T, Comp> &vec) {
    std::memcpy(mArray, (const T *)vec, sizeof(T) * Comp);
    return *this;
  }

  VectorBase &operator+() { return *this; }

  VectorBase &operator-() {
    *this *= static_cast<T>(-1);
    return *this;
  }

  VectorBase &operator+=(const VectorBase &other) {
    for (int i = 0; i < Comp; i++) {
      mArray[i] += other[i];
    }
    return *this;
  }
  VectorBase &operator-=(const VectorBase &other) {
    for (int i = 0; i < Comp; i++) {
      mArray[i] -= other[i];
    }
    return *this;
  }
  VectorBase &operator*=(const VectorBase &other) {
    for (int i = 0; i < Comp; i++) {
      mArray[i] *= other[i];
    }
    return *this;
  }
  VectorBase &operator/=(const VectorBase &other) {
    for (int i = 0; i < Comp; i++) {
      mArray[i] /= other[i];
    }
    return *this;
  }

  VectorBase &operator+=(const T &other) {
    for (int i = 0; i < Comp; i++) {
      mArray[i] += other;
    }
    return *this;
  }
  VectorBase &operator-=(const T &other) {
    for (int i = 0; i < Comp; i++) {
      mArray[i] -= other;
    }
    return *this;
  }
  VectorBase &operator*=(const T &other) {
    for (int i = 0; i < Comp; i++) {
      mArray[i] *= other;
    }
    return *this;
  }
  VectorBase &operator/=(const T &other) {
    for (int i = 0; i < Comp; i++) {
      mArray[i] /= other;
    }
    return *this;
  }

  operator T *() { return mArray; }
  operator const T *() const { return mArray; }

  T GetLength() const {
    T length = 0;
    for (int i = 0; i < Comp; i++) {
      length += mArray[i] * mArray[i];
    }
    return std::sqrt(length);
  }
  size_t GetSize() const { return mSize; }
  VectorBase GetCopy() const { return VectorBase((const T *)mArray, Comp); }
  VectorBase GetNormalized() const {
    VectorBase copy = this->GetCopy();
    copy /= this->GetLength();
    return copy;
  }

  void Normalize() { *this /= this->GetLength(); }
};

template <typename T, size_t Comp>
std::ostream &operator<<(std::ostream &stream, const VectorBase<T, Comp> &vec) {
  stream << "( ";
  for (int i = 0; i < Comp - 1; i++) {
    stream << vec[i] << ", ";
  }
  stream << vec[Comp - 1] << " )";
  return stream;
}

template <typename T, size_t Comp>
VectorBase<T, Comp> operator+(const VectorBase<T, Comp> &v1,
                              const VectorBase<T, Comp> &v2) {
  VectorBase<T, Comp> result;
  for (int i = 0; i < Comp; i++) {
    result[i] = v1[i] + v2[i];
  }
  return result;
}
template <typename T, size_t Comp>
VectorBase<T, Comp> operator-(const VectorBase<T, Comp> &v1,
                              const VectorBase<T, Comp> &v2) {
  VectorBase<T, Comp> result;
  for (int i = 0; i < Comp; i++) {
    result[i] = v1[i] - v2[i];
  }
  return result;
}
template <typename T, size_t Comp>
VectorBase<T, Comp> operator*(const VectorBase<T, Comp> &v1,
                              const VectorBase<T, Comp> &v2) {
  VectorBase<T, Comp> result;
  for (int i = 0; i < Comp; i++) {
    result[i] = v1[i] * v2[i];
  }
  return result;
}
template <typename T, size_t Comp>
VectorBase<T, Comp> operator/(const VectorBase<T, Comp> &v1,
                              const VectorBase<T, Comp> &v2) {
  VectorBase<T, Comp> result;
  for (int i = 0; i < Comp; i++) {
    result[i] = v1[i] / v2[i];
  }
  return result;
}

template <typename T, size_t Comp>
VectorBase<T, Comp> operator+(const VectorBase<T, Comp> &vec, const T &num) {
  VectorBase<T, Comp> result;
  for (int i = 0; i < Comp; i++) {
    result[i] = vec[i] + num;
  }
  return result;
}
template <typename T, size_t Comp>
VectorBase<T, Comp> operator-(const VectorBase<T, Comp> &vec, const T &num) {
  VectorBase<T, Comp> result;
  for (int i = 0; i < Comp; i++) {
    result[i] = vec[i] - num;
  }
  return result;
}
template <typename T, size_t Comp>
VectorBase<T, Comp> operator*(const VectorBase<T, Comp> &vec, const T &num) {
  VectorBase<T, Comp> result;
  for (int i = 0; i < Comp; i++) {
    result[i] = vec[i] * num;
  }
  return result;
}
template <typename T, size_t Comp>
VectorBase<T, Comp> operator/(const VectorBase<T, Comp> &vec, const T &num) {
  VectorBase<T, Comp> result;
  for (int i = 0; i < Comp; i++) {
    result[i] = vec[i] / num;
  }
  return result;
}

template <typename T, size_t Comp>
VectorBase<T, Comp> operator+(const T &num, const VectorBase<T, Comp> &vec) {
  VectorBase<T, Comp> result;
  for (int i = 0; i < Comp; i++) {
    result[i] = num + vec[i];
  }
  return result;
}
template <typename T, size_t Comp>
VectorBase<T, Comp> operator-(const T &num, const VectorBase<T, Comp> &vec) {
  VectorBase<T, Comp> result;
  for (int i = 0; i < Comp; i++) {
    result[i] = num - vec[i];
  }
  return result;
}
template <typename T, size_t Comp>
VectorBase<T, Comp> operator*(const T &num, const VectorBase<T, Comp> &vec) {
  VectorBase<T, Comp> result;
  for (int i = 0; i < Comp; i++) {
    result[i] = num * vec[i];
  }
  return result;
}
template <typename T, size_t Comp>
VectorBase<T, Comp> operator/(const T &num, const VectorBase<T, Comp> &vec) {
  VectorBase<T, Comp> result;
  for (int i = 0; i < Comp; i++) {
    result[i] = num / vec[i];
  }
  return result;
}

template <typename T, size_t Comp>
T dot(const VectorBase<T, Comp> &v1, const VectorBase<T, Comp> &v2) {
  T result = 0;
  for (int i = 0; i < 4; i++) {
    result += v1[i] * v2[i];
  }
  return result;
}

template <typename T> class vec2;
template <typename T> class vec3;
template <typename T> class vec4;

template <typename T> class vec2 : public VectorBase<T, 2> {
public:
  vec2() : VectorBase<T, 2>() { ; }
  vec2(const vec2 &data)
      : VectorBase<T, 2>((const T *)data.mArray, this->GetSize()) {
    ;
  }
  vec2(const std::vector<T> &data) : VectorBase<T, 2>(data) { ; }
  vec2(const T &c1, const T &c2) : VectorBase<T, 2>({c1, c2}) { ; }
  vec2(const vec3<T> &data)
      : VectorBase<T, 2>((const T *)data, this->GetSize()) {
    ;
  }
  vec2(const vec4<T> &data)
      : VectorBase<T, 2>((const T *)data, this->GetSize()) {
    ;
  }
};

template <typename T> class vec3 : public VectorBase<T, 3> {
public:
  vec3() : VectorBase<T, 3>() { ; }
  vec3(const vec3 &data)
      : VectorBase<T, 3>((const T *)data.mArray, this->GetSize()) {
    ;
  }
  vec3(const std::vector<T> &data) : VectorBase<T, 3>(data) { ; }
  vec3(const T &c1, const T &c2, const T &c3) : VectorBase<T, 3>({c1, c2, c3}) {
    ;
  }
  vec3(const vec2<T> &data)
      : VectorBase<T, 3>((const T *)data, data.GetSize()) {
    ;
  }
  vec3(const vec4<T> &data)
      : VectorBase<T, 3>((const T *)data, data.GetSize()) {
    ;
  }
};

template <typename T> class vec4 : public VectorBase<T, 4> {
public:
  vec4() : VectorBase<T, 4>() { ; }
  vec4(const vec4 &data)
      : VectorBase<T, 4>((const T *)data.mArray, this->GetSize()) {
    ;
  }
  vec4(const std::vector<T> &data) : VectorBase<T, 4>(data) { ; }
  vec4(const std::initializer_list<T> &data) : VectorBase<T, 4>(data) { ; }
  vec4(const T &c1, const T &c2, const T &c3, const T &c4)
      : VectorBase<T, 4>({c1, c2, c3, c4}) {
    ;
  }
  vec4(const vec2<T> &data)
      : VectorBase<T, 4>((const T *)data, data.GetSize()) {
    ;
  }
  vec4(const vec3<T> &data)
      : VectorBase<T, 4>((const T *)data, data.GetSize()) {
    ;
  }
};

template <typename T> vec3<T> cross(const vec2<T> &v1, const vec2<T> &v2) {
  vec3<T> result;
  result[2] = v1[0] * v2[1] - v1[1] * v2[0];
  return result;
}
template <typename T> vec3<T> cross(const vec3<T> &v1, const vec3<T> &v2) {
  vec3<T> result;
  result[0] = v1[1] * v2[2] - v1[2] * v2[1];
  result[1] = v1[2] * v2[0] - v1[0] * v2[2];
  result[2] = v1[0] * v2[1] - v1[1] * v2[0];
  return result;
}
} // namespace GeoMath

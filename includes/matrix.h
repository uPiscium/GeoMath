#pragma once
#include <cstring>
#include <initializer_list>
#include <iostream>
#include <ostream>
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
                        sizeof(T) * (Column < data[i].size() ? Column
                                                             : data[i].size()));
        }
    }

  public:
    MatrixBase() { this->Allocate(); }
    MatrixBase(const std::vector<std::vector<T>> &data) {
        this->Allocate();
        this->Flatten(data);
    }
    MatrixBase(const std::vector<T> &data) {
        this->Allocate();
        std::memcpy(mArray, &data[0],
                    sizeof(T) * (mSize < data.size() ? mSize : data.size()));
    }
    MatrixBase(const std::vector<VectorBase<T, Column>> &data) {
        this->Allocate();
        for (int i = 0; i < (Row < data.size() ? Row : data.size()); i++) {
            mRows[i] = data[i];
        }
    };
    MatrixBase(const T *data, const size_t &comps) {
        this->Allocate();
        std::memcpy(mArray, &data[0],
                    sizeof(T) * (mSize < comps ? mSize : comps));
    }
    MatrixBase(const MatrixBase &data) {
        this->Allocate();
        std::memcpy(mArray, &data[0], sizeof(T) * mSize);
    }
    ~MatrixBase() { this->Free(); }

    VectorBase<T, Row> &operator[](const size_t &idx) {
        if (idx < 0 || Row - 1 < idx) {
            _MATH_THROW(IndexOutOfRange("Index is out of range."));
        }
        return mRows[idx];
    }
    const VectorBase<T, Row> &operator[](const size_t &idx) const {
        if (idx < 0 || Row - 1 < idx) {
            _MATH_THROW(IndexOutOfRange("Index is out of range."));
        }
        return mRows[idx];
    }

    MatrixBase &operator=(const MatrixBase &other) {
        std::memcpy(mArray, (const T *)other, sizeof(T) * mSize);
        return *this;
    }

    MatrixBase &operator+() { return *this; }

    MatrixBase &operator-() {
        *this *= static_cast<T>(-1);
        return *this;
    }

    MatrixBase &operator+=(const MatrixBase &other) {
        for (int i = 0; i < mSize; i++) {
            mArray[i] += other[i];
        }
        return *this;
    }
    MatrixBase &operator-=(const MatrixBase &other) {
        for (int i = 0; i < mSize; i++) {
            mArray[i] -= other[i];
        }
        return *this;
    }
    MatrixBase &operator*=(const MatrixBase &other) {
        for (int i = 0; i < mSize; i++) {
            mArray[i] *= other[i];
        }
        return *this;
    }
    MatrixBase &operator/=(const MatrixBase &other) {
        for (int i = 0; i < mSize; i++) {
            mArray[i] /= other[i];
        }
        return *this;
    }

    MatrixBase &operator+=(const T &other) {
        for (int i = 0; i < mSize; i++) {
            mArray[i] += other;
        }
        return *this;
    }
    MatrixBase &operator-=(const T &other) {
        for (int i = 0; i < mSize; i++) {
            mArray[i] -= other;
        }
        return *this;
    }
    MatrixBase &operator*=(const T &other) {
        for (int i = 0; i < mSize; i++) {
            mArray[i] *= other;
        }
        return *this;
    }
    MatrixBase &operator/=(const T &other) {
        for (int i = 0; i < mSize; i++) {
            mArray[i] /= other;
        }
        return *this;
    }

    operator T *() { return mArray; }
    operator const T *() const { return mArray; }

    const size_t &GetSize() const { return mSize; }
    MatrixBase GetCopy() const {
        return MatrixBase<T, Row, Column>((const T *)mArray, mSize);
    }
    MatrixBase<T, Column, Row> GetTransposed() const {
        MatrixBase<T, Column, Row> transposed;
        for (int i = 0; i < Row; i++) {
            for (int j = 0; j < Column; j++) {
                transposed[j][i] = mArray[i * Column + j];
            }
        }
        return transposed;
    }
    MatrixBase<T, Row - 1, Column - 1> GetCofactor(const size_t &i,
                                                   const size_t &j) const {
        MatrixBase<T, Row - 1, Column - 1> result;
        int p = 0, q = 0;
        for (int l = 0; l < Row; l++) {
            if (l == i) {
                continue;
            }
            for (int m = 0; m < Column; m++) {
                if (m == j) {
                    continue;
                }
                result[p][q] = mArray[l * Row + m];
                ++q;
            }
            ++p;
        }
        return result;
    }
    // T GetDeterminant() const;
};

template <typename T, size_t Row, size_t Column>
std::ostream &operator<<(std::ostream &stream,
                         const MatrixBase<T, Row, Column> &mat) {
    stream << "( ";
    for (int i = 0; i < Row - 1; i++) {
        stream << mat[i] << std::endl;
        stream << "  ";
    }
    stream << mat[Row - 1] << " )";
    return stream;
}

template <typename T, size_t Row, size_t Column>
MatrixBase<T, Row, Column> operator+(const MatrixBase<T, Row, Column> &m1,
                                     const MatrixBase<T, Row, Column> &m2) {
    MatrixBase<T, Row, Column> result;
    for (int i = 0; i < Row; i++) {
        for (int j = 0; j < Column; j++) {
            result[i][j] = m1[i][j] + m2[i][j];
        }
    }
    return result;
}
template <typename T, size_t Row, size_t Column>
MatrixBase<T, Row, Column> operator-(const MatrixBase<T, Row, Column> &m1,
                                     const MatrixBase<T, Row, Column> &m2) {
    MatrixBase<T, Row, Column> result;
    for (int i = 0; i < Row; i++) {
        for (int j = 0; j < Column; j++) {
            result[i][j] = m1[i][j] - m2[i][j];
        }
    }
    return result;
}
template <typename T, size_t Row, size_t Column>
MatrixBase<T, Row, Column> operator*(const MatrixBase<T, Row, Column> &m1,
                                     const MatrixBase<T, Row, Column> &m2) {
    MatrixBase<T, Row, Column> result;
    for (int i = 0; i < Row; i++) {
        for (int j = 0; j < Column; j++) {
            result[i][j] = m1[i][j] * m2[i][j];
        }
    }
    return result;
}
template <typename T, size_t Row, size_t Column>
MatrixBase<T, Row, Column> operator/(const MatrixBase<T, Row, Column> &m1,
                                     const MatrixBase<T, Row, Column> &m2) {
    MatrixBase<T, Row, Column> result;
    for (int i = 0; i < Row; i++) {
        for (int j = 0; j < Column; j++) {
            result[i][j] = m1[i][j] / m2[i][j];
        }
    }
    return result;
}

template <typename T, size_t Row, size_t Column>
MatrixBase<T, Row, Column> operator+(const MatrixBase<T, Row, Column> &mat,
                                     const T &num) {
    MatrixBase<T, Row, Column> result;
    for (int i = 0; i < Row; i++) {
        for (int j = 0; j < Column; j++) {
            result[i][j] = mat[i][j] + num;
        }
    }
    return result;
}
template <typename T, size_t Row, size_t Column>
MatrixBase<T, Row, Column> operator-(const MatrixBase<T, Row, Column> &mat,
                                     const T &num) {
    MatrixBase<T, Row, Column> result;
    for (int i = 0; i < Row; i++) {
        for (int j = 0; j < Column; j++) {
            result[i][j] = mat[i][j] - num;
        }
    }
    return result;
}
template <typename T, size_t Row, size_t Column>
MatrixBase<T, Row, Column> operator*(const MatrixBase<T, Row, Column> &mat,
                                     const T &num) {
    MatrixBase<T, Row, Column> result;
    for (int i = 0; i < Row; i++) {
        for (int j = 0; j < Column; j++) {
            result[i][j] = mat[i][j] * num;
        }
    }
    return result;
}
template <typename T, size_t Row, size_t Column>
MatrixBase<T, Row, Column> operator/(const MatrixBase<T, Row, Column> &mat,
                                     const T &num) {
    MatrixBase<T, Row, Column> result;
    for (int i = 0; i < Row; i++) {
        for (int j = 0; j < Column; j++) {
            result[i][j] = mat[i][j] / num;
        }
    }
    return result;
}

template <typename T, size_t Row, size_t Column>
MatrixBase<T, Row, Column> operator+(const T &num,
                                     const MatrixBase<T, Row, Column> &mat) {
    MatrixBase<T, Row, Column> result;
    for (int i = 0; i < Row; i++) {
        for (int j = 0; j < Column; j++) {
            result[i][j] = num + mat[i][j];
        }
    }
    return result;
}
template <typename T, size_t Row, size_t Column>
MatrixBase<T, Row, Column> operator-(const T &num,
                                     const MatrixBase<T, Row, Column> &mat) {
    MatrixBase<T, Row, Column> result;
    for (int i = 0; i < Row; i++) {
        for (int j = 0; j < Column; j++) {
            result[i][j] = num - mat[i][j];
        }
    }
    return result;
}
template <typename T, size_t Row, size_t Column>
MatrixBase<T, Row, Column> operator*(const T &num,
                                     const MatrixBase<T, Row, Column> &mat) {
    MatrixBase<T, Row, Column> result;
    for (int i = 0; i < Row; i++) {
        for (int j = 0; j < Column; j++) {
            result[i][j] = num * mat[i][j];
        }
    }
    return result;
}
template <typename T, size_t Row, size_t Column>
MatrixBase<T, Row, Column> operator/(const T &num,
                                     const MatrixBase<T, Row, Column> &mat) {
    MatrixBase<T, Row, Column> result;
    for (int i = 0; i < Row; i++) {
        for (int j = 0; j < Column; j++) {
            result[i][j] = num / mat[i][j];
        }
    }
    return result;
}

template <typename T, size_t Row, size_t Share, size_t Column>
MatrixBase<T, Row, Column> dot(const MatrixBase<T, Row, Share> &m1,
                               const MatrixBase<T, Share, Column> &m2) {
    MatrixBase<T, Column, Row> tm2 = m2.GetTransposed();
    MatrixBase<T, Row, Column> result;
    for (int i = 0; i < Row; i++) {
        for (int j = 0; j < Column; j++) {
            result[i][j] = dot(m1[i], tm2[i]);
        }
    }
    return result;
}
template <typename T, size_t Row, size_t Column>
VectorBase<T, Row> dot(const MatrixBase<T, Row, Column> &mat,
                       const VectorBase<T, Column> &vec) {
    VectorBase<T, Row> result;
    for (int i = 0; i < Row; i++) {
        result[i] = dot(mat[i], vec);
    }
    return result;
}
template <typename T, size_t Row, size_t Column>
VectorBase<T, Column> dot(const VectorBase<T, Row> &vec,
                          const MatrixBase<T, Row, Column> &mat) {
    VectorBase<T, Column> result;
    for (int i = 0; i < Row; i++) {
        result[i] = dot(vec, mat[i]);
    }
    return result;
}

template <typename T, size_t Size>
T determinant(const MatrixBase<T, Size, Size> &mat) {
    T det = 0;
    for (int i = 0; i < Size; i++) {
        MatrixBase<T, Size - 1, Size - 1> &cofactor = mat.GetCofactor(i, 0);
        det += mat[i][0] * determinant(cofactor);
    }
    return det;
}
template <typename T> T determinant(const MatrixBase<T, 2, 2> &mat) {
    return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
}
template <typename T> T determinant(const MatrixBase<T, 3, 3> &mat) {
    T l1 = mat[0][0] * mat[1][1] * mat[2][2];
    T l2 = mat[0][1] * mat[1][2] * mat[2][0];
    T l3 = mat[0][2] * mat[1][0] * mat[2][1];
    T left = l1 + l2 + l3;
    T r1 = mat[0][0] * mat[1][2] * mat[2][1];
    T r2 = mat[0][1] * mat[1][0] * mat[2][2];
    T r3 = mat[0][2] * mat[1][1] * mat[2][0];
    T right = r1 + r2 + r3;
    return left - right;
}

template <typename T> class mat2;
template <typename T> class mat2x3;
template <typename T> class mat2x4;
template <typename T> class mat3x2;
template <typename T> class mat3;
template <typename T> class mat3x4;
template <typename T> class mat4x2;
template <typename T> class mat4x3;
template <typename T> class mat4;

template <typename T> class mat2 : public MatrixBase<T, 2, 2> {
  public:
    mat2() : MatrixBase<T, 2, 2>() {}
    mat2(const T *array) : MatrixBase<T, 2, 2>(array, this->GetSize()) {}
    mat2(const MatrixBase<T, 2, 2> &mat) : MatrixBase<T, 2, 2>(mat) {}
    mat2(const mat2 &mat) : MatrixBase<T, 2, 2>(mat) {}
    mat2(const vec2<T> &v1, const vec2<T> &v2)
        : MatrixBase<T, 2, 2>({v1, v2}) {}
    mat2(const T &a, const T &b, const T &c, const T &d)
        : MatrixBase<T, 2, 2>({a, b, c, d}) {}
};
template <typename T> class mat2x3 : public MatrixBase<T, 2, 3> {
  public:
    mat2x3() : MatrixBase<T, 2, 3>() {}
    mat2x3(const T *array) : MatrixBase<T, 2, 3>(array, this->GetSize()) {}
    mat2x3(const MatrixBase<T, 2, 3> &mat) : MatrixBase<T, 2, 3>(mat) {}
    mat2x3(const mat2x3 &mat) : MatrixBase<T, 2, 3>(mat) {}
    mat2x3(const vec3<T> &v1, const vec3<T> &v2)
        : MatrixBase<T, 2, 3>({v1, v2}) {}
    mat2x3(const T &a, const T &b, const T &c, const T &d, const T &e,
           const T &f)
        : MatrixBase<T, 2, 3>({a, b, c, d, e, f}) {}
};
template <typename T> class mat2x4 : public MatrixBase<T, 2, 4> {
  public:
    mat2x4() : MatrixBase<T, 2, 4>() {}
    mat2x4(const T *array) : MatrixBase<T, 2, 4>(array, this->GetSize()) {}
    mat2x4(const MatrixBase<T, 2, 4> &mat) : MatrixBase<T, 2, 4>(mat) {}
    mat2x4(const mat2x4 &mat) : MatrixBase<T, 2, 4>(mat) {}
    mat2x4(const vec4<T> &v1, const vec4<T> &v2)
        : MatrixBase<T, 2, 4>({v1, v2}) {}
    mat2x4(const T &a, const T &b, const T &c, const T &d, const T &e,
           const T &f, const T &g, const T &h)
        : MatrixBase<T, 2, 4>({a, b, c, d, e, f, g, h}) {}
};
template <typename T> class mat3x2 : public MatrixBase<T, 3, 2> {
  public:
    mat3x2() : MatrixBase<T, 3, 2>() {}
    mat3x2(const T *array) : MatrixBase<T, 3, 2>(array, this->GetSize()) {}
    mat3x2(const MatrixBase<T, 3, 2> &mat) : MatrixBase<T, 3, 2>(mat) {}
    mat3x2(const mat3x2 &mat) : MatrixBase<T, 3, 2>(mat) {}
    mat3x2(const vec2<T> &v1, const vec2<T> &v2, const vec2<T> &v3)
        : MatrixBase<T, 3, 2>({v1, v2, v3}) {}
    mat3x2(const T &a, const T &b, const T &c, const T &d, const T &e,
           const T &f)
        : MatrixBase<T, 3, 2>({a, b, c, d, e, f}) {}
};
template <typename T> class mat3 : public MatrixBase<T, 3, 3> {
  public:
    mat3() : MatrixBase<T, 3, 3>() {}
    mat3(const T *array) : MatrixBase<T, 3, 3>(array, this->GetSize()) {}
    mat3(const MatrixBase<T, 3, 3> &mat) : MatrixBase<T, 3, 3>(mat) {}
    mat3(const mat3 &mat) : MatrixBase<T, 3, 3>(mat) {}
    mat3(const vec3<T> &v1, const vec3<T> &v2, const vec3<T> &v3)
        : MatrixBase<T, 3, 3>({v1, v2, v3}) {}
    mat3(const T &a, const T &b, const T &c, const T &d, const T &e, const T &f,
         const T &g, const T &h, const T &i)
        : MatrixBase<T, 3, 3>({a, b, c, d, e, f, g, h, i}) {}
};
template <typename T> class mat3x4 : public MatrixBase<T, 3, 4> {
  public:
    mat3x4() : MatrixBase<T, 3, 4>() {}
    mat3x4(const T *array) : MatrixBase<T, 3, 4>(array, this->GetSize()) {}
    mat3x4(const MatrixBase<T, 3, 4> &mat) : MatrixBase<T, 3, 4>(mat) {}
    mat3x4(const mat3x4 &mat) : MatrixBase<T, 3, 4>(mat) {}
    mat3x4(const vec4<T> &v1, const vec4<T> &v2, const vec4<T> &v3)
        : MatrixBase<T, 3, 4>({v1, v2, v3}) {}
    mat3x4(const T &a, const T &b, const T &c, const T &d, const T &e,
           const T &f, const T &g, const T &h, const T &i, const T &j,
           const T &k, const T &l)
        : MatrixBase<T, 3, 4>({a, b, c, d, e, f, g, h, i, j, k, l}) {}
};
template <typename T> class mat4x2 : public MatrixBase<T, 4, 2> {
  public:
    mat4x2() : MatrixBase<T, 4, 2>() {}
    mat4x2(const T *array) : MatrixBase<T, 4, 2>(array, this->GetSize()) {}
    mat4x2(const MatrixBase<T, 4, 2> &mat) : MatrixBase<T, 4, 2>(mat) {}
    mat4x2(const mat4x2 &mat) : MatrixBase<T, 4, 2>(mat) {}
    mat4x2(const vec2<T> &v1, const vec2<T> &v2, const vec2<T> &v3,
           const vec2<T> &v4)
        : MatrixBase<T, 4, 2>({v1, v2, v3, v4}) {}
    mat4x2(const T &a, const T &b, const T &c, const T &d, const T &e,
           const T &f, const T &g, const T &h)
        : MatrixBase<T, 4, 2>({a, b, c, d, e, f, g, h}) {}
};
template <typename T> class mat4x3 : public MatrixBase<T, 4, 3> {
  public:
    mat4x3() : MatrixBase<T, 4, 3>() {}
    mat4x3(const T *array) : MatrixBase<T, 4, 3>(array, this->GetSize()) {}
    mat4x3(const MatrixBase<T, 4, 3> &mat) : MatrixBase<T, 4, 3>(mat) {}
    mat4x3(const mat4x3 &mat) : MatrixBase<T, 4, 3>(mat) {}
    mat4x3(const vec3<T> &v1, const vec3<T> &v2, const vec3<T> &v3,
           const vec3<T> &v4)
        : MatrixBase<T, 4, 3>({v1, v2, v3, v4}) {}
    mat4x3(const T &a, const T &b, const T &c, const T &d, const T &e,
           const T &f, const T &g, const T &h, const T &i, const T &j,
           const T &k, const T &l)
        : MatrixBase<T, 4, 3>({a, b, c, d, e, f, g, h, i, j, k, l}) {}
};
template <typename T> class mat4 : public MatrixBase<T, 4, 4> {
  public:
    mat4() : MatrixBase<T, 4, 4>() {}
    mat4(const T *array) : MatrixBase<T, 4, 4>(array, this->GetSize()) {}
    mat4(const MatrixBase<T, 4, 4> &mat) : MatrixBase<T, 4, 4>(mat) {}
    mat4(const mat4 &mat) : MatrixBase<T, 4, 4>(mat) {}
    mat4(const vec4<T> &v1, const vec4<T> &v2, const vec4<T> &v3,
         const vec4<T> &v4)
        : MatrixBase<T, 4, 4>({v1, v2, v3, v4}) {}
    mat4(const T &a, const T &b, const T &c, const T &d, const T &e, const T &f,
         const T &g, const T &h, const T &i, const T &j, const T &k, const T &l,
         const T &m, const T &n, const T &o, const T &p)
        : MatrixBase<T, 4, 4>(
              {a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p}) {}
};
} // namespace GeoMath

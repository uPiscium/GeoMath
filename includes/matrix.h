#include <iostream>
#include <initializer_list>
#include <vector>

#include "exceptions.h"
namespace GeoMath
{
	template <typename T, size_t Row, size_t Column>
	class Matrix
	{
	private:
		using value_type = T;
		using size_type = size_t;

	private:
		size_t mSize = Row * Column;
		T *mArray = nullptr;
		bool mSelfAlloc = false;

	private:
		void Allocate()
		{
			mArray = new T[mSize];
			for (int i = 0; i < mSize; i++)
			{
				mArray[i] = static_cast<T>(0);
			}
			mSelfAlloc = true;
		}
		void Free()
		{
			if (mSelfAlloc)
			{
				delete[] mArray;
			}
		}

	public:
		Matrix() { this->Allocate(); }
		Matrix(const std::vector<T> &data)
		{
			this->Allocate();
			size_t copysize = sizeof(T) * (mSize <= data.size() ? mSize : data.size());
			std::memcpy(mArray, &data[0], copysize);
		}
		Matrix(const std::initializer_list<T> &data)
		{
			this->Allocate();
			size_t copysize = sizeof(T) * (mSize <= data.size() ? mSize : data.size());
			std::memcpy(mArray, &data[0], copysize);
		}
		Matrix(const std::vector<Matrix<T, 1, Column>> &data)
		{
			for (int i = 0; i < data.size(); i++)
			{
				std::memcpy(&mArray[i * Column])
			}
		}
		Matrix(const std::vector<Matrix<T, Row, 1>> &data);
		Matrix(const Matrix &data);

		T operator[](const size_t &idx)
		{
			if (idx < 0 || mSize - 1 < idx)
			{
				_MATH_THROW(IndexOutOfRange("Index is out of range."));
			}
			return mArray[idx];
		}
	};
}

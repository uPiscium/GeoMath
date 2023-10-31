#include <iostream>
#include <initializer_list>
#include <vector>

#include "exceptions.h"

namespace GeoMath
{
    template <typename T, size_t Comp>
    class Vector
    {
    private:
        using value_type = T;

    private:
        size_t mSize = Comp;
        T *mArray = nullptr;
        bool mSelfAlloc = false;

    private:
        void Allocate()
        {
            mArray = new T[mSize];
            mSelfAlloc = true;
        }

    public:
        Vector();
        Vector(const std::vector<T>& data);
        Vector(const std::initializer_list<T>& data);
        template <typename std::enable_if<Comp == 1, nullptr> = nullptr>
        Vector(T c1);
        template <typename std::enable_if<Comp == 2, nullptr> = nullptr>
        Vector(T c1, T c2);
        template <typename std::enable_if<Comp == 3, nullptr> = nullptr>
        Vector(T c1, T c2, T c3);
        template <typename std::enable_if<Comp == 4, nullptr> = nullptr>
        Vector(T c1, T c2, T c3, T c4);
    };
}
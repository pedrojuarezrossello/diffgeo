#ifndef MATHS_H
#define MATHS_H
#include <cfloat>
#include <cmath>

constexpr auto M_PI = 3.14159265358979323846;

template<typename T,typename U,
         std::enable_if_t<std::is_floating_point_v<T> && std::is_floating_point_v<U>, std::nullptr_t> = nullptr>
bool almostEqualRelativeAndAbs(T A, U B,
    double maxDiff, float maxRelDiff = FLT_EPSILON)
{
    // Check if the numbers are really close -- needed
    // when comparing numbers near zero.
    auto diff = fabs(A - B);
    if (diff <= maxDiff)
        return true;

    A = fabs(A);
    B = fabs(B);

    auto largest = (B > A) ? B : A;

    if (diff <= largest * maxRelDiff)
        return true;
    return false;
}

template<typename T,
    std::enable_if_t<std::is_floating_point_v<T>, std::nullptr_t> = nullptr>
bool isZero(T a)
{
    return almostEqualRelativeAndAbs(a, 0.0, 1.0e-6);
}

template<typename T,
    std::enable_if_t<std::is_floating_point_v<T>, std::nullptr_t> = nullptr>
auto signOf(T num)
{
    return (num >= 0) ? "+" : "";
}
#endif //!MATHS_H
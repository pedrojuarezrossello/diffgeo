#ifndef LINE_H
#define LINE_H
#include "../vector/vector.h"

template<typename T,
    std::enable_if_t<std::is_floating_point_v<T>, std::nullptr_t> = nullptr>
class Line
{
	Vector<T> direction_;
    Vector<T> point_;
public:
    Line( Vector<T>& point, Vector<T>& direction)
						: direction_(direction), point_(point) {}
    const Vector<T>& getDirection() const
    {
        return direction_;
    }
	Vector<T>& getDirection()
    {
        return const_cast<int&>(const_cast<const Line*>(this)->getDirection());
    }
    template<typename U,
        std::enable_if_t<std::is_floating_point_v<T>, std::nullptr_t> = nullptr >
    auto at(U x)
    {
        auto newPoint = point_ + x * direction_;
        return newPoint; //(N)RVO
    }
    template<typename U,
        std::enable_if_t<std::is_floating_point_v<U>, std::nullptr_t> = nullptr>
    bool isInLine(const Vector<U>& vec)
    {
        return areParallel(vec - point_, direction_);
    }
};


#endif //!LINE_H

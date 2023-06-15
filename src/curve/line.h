#ifndef LINE_H
#define LINE_H
#include "../vector/vector.h"
#include "../utils/operations.h"

template<typename T, typename U,
    std::enable_if_t<std::is_floating_point_v<T> && std::is_floating_point_v<U>, std::nullptr_t> = nullptr>
class Line
{
	Vector<T> direction_;
    Vector<U> point_;

public:

    //constructors
    Line(const Vector<U>& point, const Vector<T>& direction)
						: direction_(direction), point_(point)
    {
        direction_.normalise();
    }

    Line(const Vector<T>& point_a, const Vector<T>& point_b, bool pos_direction)
						: direction_(pos_direction ? (point_b-point_a):(point_a-point_b)), point_(point_a)
    {
        direction_.normalise();
    }

    //getters
    const Vector<T>& getDirection() const
    {
        return direction_;
    }
	Vector<T>& getDirection()
    {
        return const_cast<int&>(const_cast<const Line*>(this)->getDirection());
    }

    //print util
    friend ostream& operator<<(ostream& os, Line const& line)
    {
        os << line.point_ << " + t" << line.direction_;
        return os;
    }

    //fetch point at x
    template<typename V,
        std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
    auto at(V x)
    {
        auto newPoint = point_ + x * direction_;
        return newPoint; //(N)RVO
    }

    //check for inclusion of a point 
    template<typename V,
        std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr>
    bool isInLine(const Vector<V>& vec)
    {
        return areParallel(vec - point_, direction_);
    }

    //parallel lines
    template<typename V,typename W,
        std::enable_if_t<std::is_floating_point_v<V> && std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
    friend bool areParallel(const Line& line1, const Line<V,W>& line2)
    {
        return areParallel(line1.direction_, line2.direction_);
    }

    //distance between two lines
    template<typename V, typename W,
        std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
    friend double distance(const Line& line1, const Line<V,W>& line2)
    {
        if (!areParallel(line1.direction_, line2.direction_)) 
        {
            auto normalVector(cross_product(line1.direction_, line2.direction_));
            normalVector.normalise();
            return fabs(normalVector * (line1.point_ - line2.point_));
        } 

    	auto normalVector = cross_product(line1.point_ - line2.point_, line1.direction_);
    	return normalVector.norm();
    }

    //closest point on line from other point
    template<typename Point,
        typename = std::enable_if_t<is_vector_or_expression_t<Point>>>
    auto closestPoint(const Point& point)
    {
        const double projection = (point - point_) * direction_;
        auto closestPoint = at(projection);
        return closestPoint;
    }

    //distance between point and line
    template<typename U,
        std::enable_if_t<std::is_floating_point_v<U>, std::nullptr_t> = nullptr>
    double distanceFrom(const Vector<U>& point)
    {
        auto closest(closestPoint(point));
        return distance(closest, point);
    }

    template<typename CallableObject, typename... Args>
	double distanceFrom(const VectorExpression<CallableObject, Args...>& point)
    {
        const Vector<double> pointAsVec = point;
        return distanceFrom(pointAsVec);
    }

    //perpendicular line through point P
    template<typename V,
        std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr>
    Line<V,T> perpendicularFrom(const Vector<V>& point)
    {
        Line<V,T> perpendicular(point, direction_.perpendicular());
        return perpendicular; //(N)RVO
    }

    //parallel line through point P
    template<typename V,
        std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr>
    Line<V, T> parallelFrom(const Vector<V>& point)
    {
        Line<V, T> parallel(point, direction_);
        return parallel; //(N)RVO
    }
};

	

#endif //!LINE_H

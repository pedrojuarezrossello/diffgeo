#ifndef DG_LINE_H
#define DG_LINE_H
#include "../vector/vector.h"
#include "../vector/operations.h"

namespace dg {

    namespace curve {

        template<typename T, typename U,
            std::enable_if_t<std::is_floating_point_v<T>&& std::is_floating_point_v<U>, std::nullptr_t> = nullptr>
        class Line
        {
            dg::vector::Vector<T> direction_;
            dg::vector::Vector<U> point_;

        public:

            //Constructors 

            Line() = default;

            Line(const dg::vector::Vector<U>& point, const dg::vector::Vector<T>& direction)
                : direction_(direction), point_(point)
            {
                direction_.normalise();
            }

            template<typename CallableObject, typename ... Args>
            Line(const dg::vector::Vector<U>& point, const dg::vector::VectorExpression<CallableObject, Args...>& direction)
                : direction_(direction), point_(point)
            {
                direction_.normalise();
            }

            template<typename CallableObject, typename ... Args>
            Line(const dg::vector::VectorExpression<CallableObject, Args...>& point, const dg::vector::Vector<U>& direction)
                : direction_(direction), point_(point)
            {
                direction_.normalise();
            }

            template<typename CallableObject1, typename ... Args1,
                typename CallableObject2, typename ...Args2>
            Line(const dg::vector::VectorExpression<CallableObject1, Args1...>& point, const dg::vector::VectorExpression<CallableObject2, Args2...>& direction)
                : direction_(direction), point_(point)
            {
                direction_.normalise();
            }

            Line(const dg::vector::Vector<T>& point_a, const dg::vector::Vector<T>& point_b, bool pos_direction)
                : direction_(pos_direction ? (point_b - point_a) : (point_a - point_b)), point_(point_a)
            {
                direction_.normalise();
            }

            //Equality

            bool operator==(const Line& rhs)
            {
                return direction_ == rhs.direction_ && point_ == rhs.point_;
            }

            //Getter
            const dg::vector::Vector<T>& getDirection() const
            {
                return direction_;
            }

            //Print util
            friend std::ostream& operator<<(std::ostream& os, Line const& line)
            {
                os << line.point_ << " + t" << line.direction_;
                return os;
            }

            //Fetch point at x
            template<typename V,
                std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
            auto at(V x) const
            {
                dg::vector::Vector<T> newPoint(point_ + x * direction_);
                return newPoint; //(N)RVO
            }

            //Check for inclusion of a point 
            template<typename V,
                std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr>
            bool isInLine(const dg::vector::Vector<V>& vec)
            {
                auto pointToPoint(dg::vector::Vector<T>(vec - point_));
                return areParallel(direction_, pointToPoint);
            }

            //Parallel lines
            template<typename V, typename W,
                std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
            friend bool areParallel(const Line& line1, const Line<V, W>& line2)
            {
                return areParallel(line1.direction_, line2.direction_);
            }

            //Parallel line through point P
            template<typename V,
                std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr>
            Line<V, T> parallelFrom(const dg::vector::Vector<V>& point)
            {
                Line<V, T> parallel(point, direction_);
                return parallel; //(N)RVO
            }

            //Intersection
            template<typename V, typename W,
                std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
            bool intersects(const Line<V, W>& line)
            {
                if (areParallel(*this, line) && !*this == line) return false;
                //two lines contained in the same plane
                return dg::math::isZero(dg::vector::cross_product(direction_, line.direction_) * (point_ - line.point_));

            }

            //Distance between two lines
            template<typename V, typename W,
                std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
            friend double distance(const Line& line1, const Line<V, W>& line2)
            {
                if (!areParallel(line1, line2))
                {
                    auto normalVector(dg::vector::cross_product(line1.direction_, line2.direction_));
                    normalVector.normalise();
                    return fabs(normalVector * (line1.point_ - line2.point_));
                }

                auto normalVector(dg::vector::cross_product(line1.point_ - line2.point_, line1.direction_));
                return normalVector.norm();
            }

            //Closest point on line from other point
            template<typename Point,
                typename = std::enable_if_t<dg::vector::is_vector_or_expression_t<Point>>>
            auto closestPoint(const Point& point)
            {
                const double projection = (point - point_) * direction_;
                auto closestPoint = at(projection);
                return closestPoint;
            }

            //Distance between point and line
            template<typename V,
                std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr>
            double distanceFrom(const dg::vector::Vector<V>& point)
            {
                auto closest(closestPoint(point));
                return distance(closest, point);
            }

            template<typename CallableObject, typename... Args>
            double distanceFrom(const dg::vector::VectorExpression<CallableObject, Args...>& point)
            {
                const dg::vector::Vector<double> pointAsVec(point);
                return distanceFrom(pointAsVec);
            }

            //Perpendicular line through point P
            template<typename V,
                std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr>
            Line<V, T> perpendicularFrom(const dg::vector::Vector<V>& point)
            {
                Line<V, T> perpendicular(point, direction_.perpendicular());
                return perpendicular; //(N)RVO
            }

           
            //Angle between two lines
            template<typename V, typename W,
                std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
            friend double angle(const Line& line1, const Line<V, W>& line2)
            {
                return angle(line1.direction_, line2.direction_);
            }

           
        };

    } //namespace curve
    
} //namespace dg

#endif //!DG_LINE_H

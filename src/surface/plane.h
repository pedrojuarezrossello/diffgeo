#ifndef PLANE_H
#define PLANE_H
#include <cstddef>
#include "../vector/vector.h"
#include "../utils/operations.h"
#include "../curve/line.h"
#include "../utils/maths.h"

namespace dg {

    namespace surf {

        template<typename T, typename U,
            std::enable_if_t<std::is_floating_point_v<T>&& std::is_floating_point_v<U>, std::nullptr_t> = nullptr>
        class Plane
        {
            dg::vector::Vector<T> normal_;
            dg::vector::Vector<U> point_;

        public:

            //constructors & equality

            Plane() = default;

            Plane(const dg::vector::Vector<T>& normal, const dg::vector::Vector<U>& point)
                : normal_(normal), point_(point)
            {
                normal_.normalise();
            }

            Plane(const dg::vector::Vector<T>& vec1, const dg::vector::Vector<T>& vec2, const dg::vector::Vector<U>& point)
                : normal_(dg::vector::cross_product(vec1, vec2)), point_(point)
            {
                normal_.normalise();
            }

            template<typename V,
                std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr>
            Plane(const dg::curve::Line<T, V>& line, const dg::vector::Vector<U>& point)
                : normal_(dg::vector::cross_product(line.getDirection(), line.at(0.0) - point)), point_(point)
            {
                normal_.normalise();
            }

            template<typename V, typename W,
                std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
            Plane(const dg::curve::Line<T, V>& line1, const dg::curve::Line<U, W>& line2)
                : normal_(dg::vector::Vector<T>()), point_(dg::vector::Vector<U>())
            {
                using namespace dg::curve;
                if (line1.intersects(line2))
                {
                    normal_ = dg::vector::cross_product(line1.getDirection(), line2.getDirection);  normal_.normalise();
                    point_ = line1.at(0.0);
                }
                else if (areParallel(line1, line2))
                {
                    normal_ = dg::vector::cross_product(line1.getDirection(), line1.at(0.0) - line2.at(0.0)); normal_.normalise();
                    point_ = line1.at(0.0);
                }
                else { throw "Lines must intersect or be parallel"; }
            }

            template <typename CallableObject, typename... Args>
            Plane(const dg::vector::VectorExpression<CallableObject, Args...>& vector_expression, const dg::vector::Vector<U>& point)
                : normal_(vector_expression), point_(point)
            {
                normal_.normalise();
            }

            template <typename CallableObject, typename... Args>
            Plane(const dg::vector::Vector<T>& normal, const dg::vector::VectorExpression<CallableObject, Args...>& point)
                : normal_(normal), point_(point)
            {
                normal_.normalise();
            }

            template <typename CallableObject1, typename... Args1, typename CallableObject2, typename... Args2>
            Plane(const dg::vector::VectorExpression<CallableObject1, Args1...>& vector_expression, const dg::vector::VectorExpression<CallableObject2, Args2...>& point)
                : normal_(vector_expression), point_(point)
            {
                normal_.normalise();
            }

            bool operator==(const Plane& other) const
            {
                return normal_ == other.normal_ && point_ == other.point_;
            }

            //getters
            const dg::vector::Vector<T>& getNormal() const { return normal_; }

            const dg::vector::Vector<T>& getPoint() const { return point_; }

            //print util
            friend std::ostream& operator<<(std::ostream& os, Plane const& plane)
            {
                auto sign{ "" };

                const T scalar{ plane.normal_ * plane.point_ };
                if (!dg::math::isZero(plane.normal_[0]))
                {
                    os << plane.normal_[0] << "x";
                }
                if (!dg::math::isZero(plane.normal_[1]))
                {
                    sign = dg::math::signOf(plane.normal_[1]);
                    os << sign << plane.normal_[1] << "y";
                }
                if (!dg::math::isZero(plane.normal_[2]))
                {
                    sign = dg::math::signOf(plane.normal_[2]);
                    os << sign << plane.normal_[2] << "z";
                }
                sign = dg::math::signOf(scalar);
                os << sign << scalar << "=0";
                return os;
            }

            //is a point on the plane?
            template<typename V,
                std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr>
            bool isOnPlane(const dg::vector::Vector<V>& point)
            {
                return dg::math::isZero(normal_ * (point - point_));
            }

            //is a line on the plane?
            template<typename V, typename W,
                std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
            bool isOnPlane(const dg::curve::Line<V, W>& line)
            {
                return isOnPlane(line.at(0.0)) && isOnPlane(line.at(1.0));
            }

            //projections
            template<typename V,
                std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr>
            auto projection(const dg::vector::Vector<V>& point)
            {
                dg::vector::Vector<V> projOnNormal = dg::vector::Vector((point * normal_) * normal_);
                dg::vector::Vector<V> projOnPlane = point - projOnNormal;
                return projOnPlane;
            }

            template<typename V, typename W,
                std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
            auto projection(const dg::curve::Line<V, W>& line)
            {
                auto vec = dg::vector::Vector(this->projection(line.getDirection()));
                auto point = this->intersectsAt(line);
                auto projLine = dg::curve::Line(point, vec);
                return projLine;
            }

            //distance
            template<typename VectorType,
                std::enable_if_t<dg::vector::is_vector_or_expression_t<VectorType>, std::nullptr_t> = nullptr>
            double distanceFrom(const VectorType& point) const
            {
                return fabs((point-point_) * normal_);
            }

            template<typename VectorType,
                std::enable_if_t<dg::vector::is_vector_or_expression_t<VectorType>, std::nullptr_t> = nullptr>
            double signedDistanceFrom(const VectorType& point) const
            {
                return (point - point_) * normal_;
            }

            template<typename V, typename W,
                std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
            double distance(const dg::curve::Line<V, W>& line) 
            {
                return areParallel(*this, line) ? this->distanceFrom(line.at(0)) : 0.0;
            }

            template<typename V, typename W,
                std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
                double distance(const Plane<V, W>& plane)
            {
                return areParallel(*this, plane) ? this->distanceFrom(plane.getPoint()) : 0.0;
            }
    
            //parallel
            template<typename V, typename W,
                std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
            friend bool areParallel(const Plane& plane1, const Plane<V, W>& plane2)
            {
                return plane1.normal_ == plane2.normal_;
            }

            template<typename V, typename W,
                std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
            friend bool areParallel(Plane& plane1, const Plane<V, W>& plane2)
            {
                return plane1.normal_ == plane2.normal_;
            }

            template<typename V, typename W,
                std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
            friend bool areParallel(const Plane& plane, const dg::curve::Line<V, W>& line)
            {
                return plane.normal_.isPerpendicular(line.getDirection());
            }

            template<typename V, typename W,
                std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
            friend bool areParallel(Plane& plane, const dg::curve::Line<V, W>& line)
            {
                return plane.normal_.isPerpendicular(line.getDirection());
            }

            //intersection
            template<typename V, typename W,
                std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
            bool intersects(const dg::curve::Line<V, W>& line)
            {
                return !areParallel(*this, line);
            }

            template<typename V, typename W,
                std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
            bool intersects(const Plane<V, W>& plane)
            {
                return !areParallel(*this, plane);
            }

            template<typename V, typename W,
                std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
            auto intersectsAt(const dg::curve::Line<V, W>& line)
            {
                if (!this->intersects(line)) throw std::exception("Lines do not intersect!");

                const auto parameter{ ((point_ - line.at(0.0)) * normal_) / (line.getDirection() * normal_) };
                auto pointIntersection = dg::vector::Vector(line.at(parameter));
                return pointIntersection;

            }

            template<typename V, typename W,
                std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
            auto intersectsAlong(const Plane<V, W>& plane)
            {
                if (!this->intersects(plane)) throw std::exception("Planes don't intersect!");

                const auto planeConst1{ normal_ * point_ };
                const auto planeConst2{ plane.normal_ * plane.point_ };
                const auto dottedNormals{ normal_ * plane.normal_ };

                const auto coeff1{ (planeConst1 - planeConst2 * dottedNormals) / (1 - dottedNormals * dottedNormals) };
                const auto coeff2{ (planeConst2 - planeConst1 * dottedNormals) / (1 - dottedNormals * dottedNormals) };

                auto lineIntersection(dg::curve::Line<T, U>(coeff1 * normal_ + coeff2 * plane.normal_, dg::vector::cross_product(normal_, plane.normal_)));
                return lineIntersection;

            }

        };

    } //namespace surf

} //namespace dg

#endif //!PLANE_H
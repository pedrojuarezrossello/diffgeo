#ifndef CIRCUMFERENCE_H
#define CIRCUMFERENCE_H
#include "../surface/plane.h"

namespace dg {

    namespace curve {

        template<typename T, typename U, typename V,
            std::enable_if_t<std::is_floating_point_v<T>&& std::is_floating_point_v<U>&& std::is_floating_point_v<V>, std::nullptr_t> = nullptr>
        class Circumference
        {
            T radius_;
            dg::vector::Vector<T> centre_;
            dg::surf::Plane<U, V> plane_;

        public:

            //constructor & equality
            Circumference() = default;

            Circumference(const T& radius, const dg::vector::Vector<T>& centre, const dg::surf::Plane<U, V>& plane)
                : radius_(radius), centre_(centre), plane_(plane) {}

            template<typename CallableObject, typename ... Args>
            Circumference(const T& radius, const dg::vector::VectorExpression<CallableObject, Args...> centre, const dg::surf::Plane<U, V>& plane)
                : radius_(radius), centre_(centre), plane_(plane) {}

            bool operator==(const Circumference& other) const
            {
                return dg::math::almostEqualRelativeAndAbs(radius_, other.radius_, 1.0e-6) && centre_ == other.centre_ && plane_ = other.plane_;
            }

            //getters 
            const T& getRadius() const { return radius_; }
            const dg::vector::Vector<T>& getCentre() const { return centre_; }
            const dg::plane::Plane<U, V>& getPlane() const { return plane_; }

            //point
            dg::vector::Vector<T> at(T theta)
            {
                const Vector<U> normal(plane_.getNormal());
                const Vector<U> plane_vector1(normal.perpendicular());
                Vector<U> plane_vector2(cross_product(normal, plane_vector1));
                plane_vector2.normalise();

                Vector<T> point(centre_ + radius_ * sin(theta) * plane_vector1 + radius_ * cos(theta) * plane_vector2);
                return point;
            }

            //length
            T length() { return 2 * dg::math::PI * radius_; }
        };

    } //namespace curve

} //namespace dg







#endif //!CIRCUMFERENCE_H
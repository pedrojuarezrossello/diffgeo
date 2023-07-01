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
            dg::surf::Plane<U, V> plane_;

        public:

            //constructor & equality
            Circumference() = default;

            Circumference(const T& radius, const dg::surf::Plane<U, V>& plane)
                : radius_(radius), plane_(plane) {}

            bool operator==(const Circumference& other) const
            {
                return dg::math::almostEqualRelativeAndAbs(radius_, other.radius_, 1.0e-6) && plane_ = other.plane_;
            }

            //getters 
            T getRadius() const { return radius_; }
            T getDiameter() const { return 2 * radius_; }
            const dg::vector::Vector<V>& getCentre() const { return plane_.getPoint(); }
            const dg::surf::Plane<U, V>& getPlane() const { return plane_; }
            
            //point
            dg::vector::Vector<T> at(T theta)
            {
                const dg::vector::Vector<U> normal(plane_.getNormal());
                const  dg::vector::Vector<U> plane_vector1(normal.perpendicular());
                dg::vector::Vector<U> plane_vector2(dg::vector::cross_product(normal, plane_vector1));
                plane_vector2.normalise();

                dg::vector::Vector<T> point(plane_.getPoint() + radius_ * sin(theta) * plane_vector1 + radius_ * cos(theta) * plane_vector2);
                return point;
            }

            //length
            T length() { return 2 * dg::math::PI * radius_; }

            //curvature
            T curvature() { return 1.0 / radius_; }

            //torsion
            T torsion() { return 0.0; }

        };

    } //namespace curve

} //namespace dg

#endif //!CIRCUMFERENCE_H
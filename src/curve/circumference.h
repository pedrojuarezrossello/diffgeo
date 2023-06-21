#ifndef CIRCUMFERENCE_H
#define CIRCUMFERENCE_H
#include "../surface/plane.h"

template<typename T, typename U, typename V,
    std::enable_if_t<std::is_floating_point_v<T>&& std::is_floating_point_v<U> && std::is_floating_point_v<V>, std::nullptr_t> = nullptr>
class Circumference 
{
     T radius_;
     Vector<T> centre_;
     Plane<U, V> plane_;

public:

    //constructor & equality
    Circumference() = default;

    Circumference(const T& radius, const Vector<T>& centre, const Plane<U, V>& plane)
        : radius_(radius), centre_(centre), plane_(plane) {}

    template<typename CallableObject, typename ... Args>
    Circumference(const T& radius, const VectorExpression<CallableObject, Args...> centre, const Plane<U, V>& plane) 
        : radius_(radius), centre_(centre), plane_(plane) {}
    
    bool operator==(const Circumference& other) const
    {
        return almostEqualRelativeAndAbs(radius_, other.radius_, 1.0e-6) && centre_ == other.centre_ && plane_ = other.plane_;
    }

    //getters 
    const T& getRadius() const { return radius_; }
    const Vector<T>& getCentre() const { return centre_; }
    const Plane<U, V>& getPlane() const {return plane_; }

    //
    Vector<T> at(T theta) 
    {
        const Vector<U> normal( plane_.getNormal() );
        const Vector<U> plane_vector1(normal.perpendicular());
        Vector<U> plane_vector2(cross_product(normal, plane_vector1));
        plane_vector2.normalise();

        Vector<T> point(centre_ + radius_*sin(theta) * plane_vector1 + radius_*cos(theta) * plane_vector2);
        return point;
    }

};








#endif //!CIRCUMFERENCE_H
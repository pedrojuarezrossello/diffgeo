#ifndef SPHERE_H
#define SPHERE_H
#include <type_traits>
#include "../vector/vector.h"


template<typename T,
         std::enable_if_t<std::is_floating_point_v<T>, std::nullptr_t> = nullptr>
class Sphere
{
	Vector<T> centre_;
	T radius_;

public:

	//constructors & equality
	Sphere(const Vector<T>& centre, T radius)
		: centre_(centre), radius_(radius) {}

	Sphere() = default;

	Sphere(T a, T b, T c, T radius)
		: centre_(a, b, c), radius_(radius){}

	template<typename CallableObject, typename ... Args>
	Sphere(const VectorExpression<CallableObject, Args...>& rhs, T radius) : centre_(rhs), radius_(radius) {}

	bool operator==(const Sphere& rhs)
	{
		return (centre_ == rhs.centre_) && almostEqualRelativeAndAbs(radius_, rhs.radius_, 1.0e-6);
	}

	//print util
	friend ostream& operator<<(ostream& os, Sphere const& sphere)
	{
		os << "||x-" << sphere.centre_ << "|| = " << sphere.radius_;
		return os;
	}

	//check if a point is on the sphere
	template<typename V,
		std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr>
	bool isOnSphere(const Vector<V>& point)
	{
		return almostEqualRelativeAndAbs(radius_, distance(centre_, point));
	}

	//surface area
	double area()
	{
		return 4 * M_PI * radius_ * radius_;
	}
};
#endif //!SPHERE_H


#ifndef SPHERE_H
#define SPHERE_H
#include <type_traits>
#include "../vector/vector.h"

namespace dg {

	namespace surf {

		template<typename T,
			std::enable_if_t<std::is_floating_point_v<T>, std::nullptr_t> = nullptr>
		class Sphere
		{
			dg::vector::Vector<T> centre_;
			T radius_;

		public:

			//constructors & equality
			Sphere(const dg::vector::Vector<T>& centre, T radius)
				: centre_(centre), radius_(radius) {}

			Sphere() = default;

			Sphere(T a, T b, T c, T radius)
				: centre_(a, b, c), radius_(radius) {}

			template<typename CallableObject, typename ... Args>
			Sphere(const dg::vector::VectorExpression<CallableObject, Args...>& rhs, T radius) : centre_(rhs), radius_(radius) {}

			bool operator==(const Sphere& rhs)
			{
				return (centre_ == rhs.centre_) && dg::math::almostEqualRelativeAndAbs(radius_, rhs.radius_, 1.0e-6);
			}

			//print util
			friend std::ostream& operator<<(std::ostream& os, Sphere const& sphere)
			{
				os << "||x-" << sphere.centre_ << "|| = " << sphere.radius_;
				return os;
			}

			//theta varies from 0 to pi and phi varies from 0 to 2pi
			dg::vector::Vector<T> at(T theta, T phi)
			{
				std::vector::Vector<T> pointOnSphere(centre_[0] + radius_ * sin(theta) * cos(phi),
					centre_[1] + radius_ * sin(theta) * sin(phi),
					centre_[2] + radius_ * cos(theta));

				return pointOnSphere;

			}

			//check if a point is on the sphere
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr>
			bool isOnSphere(const dg::vector::Vector<V>& point)
			{
				return dg::math::almostEqualRelativeAndAbs(radius_, distance(centre_, point));
			}

			//surface area
			double area()
			{
				return 4 * dg::math::PI * radius_ * radius_;
			}
		};

	} //namespace surf

} //namespace dg

#endif //!SPHERE_H


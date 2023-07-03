#ifndef SPHERE_H
#define SPHERE_H
#include <type_traits>
#include "plane.h"
#include "../curve/circumference.h"
#include <optional>

namespace dg {

	namespace surf {

		template<typename T,
			std::enable_if_t<std::is_floating_point_v<T>, std::nullptr_t> = nullptr>
		class Sphere
		{
			dg::vector::Vector<T> centre_;
			T radius_;

			//intersection and tangent helper methods
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr>
			bool isInternalTangent_(const Sphere<V>& sphere) const 
			{
				auto dist{ distance(centre_, sphere.centre_) };
				return dg::math::almostEqualRelativeAndAbs(dist + std::min(radius_, sphere.radius_), std::max(radius_, sphere.radius_), 1.0e-6);
			}

			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr>
			bool isExternalTangent_(const Sphere<V>& sphere) const 
			{
				auto dist{ distance(centre_, sphere.centre_) };
				return dg::math::almostEqualRelativeAndAbs(dist, radius_ + sphere.radius_, 1.0e-6);
			}

			template<typename V,
			std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			Plane<V, V> intersectionPlane_(const Sphere<V>& sphere) const
			{
				auto d{ distance(centre_, sphere.centre_) };
				auto h{ 0.5 + (radius_ * radius_ - sphere.radius_ * sphere.radius_) / (2.0 * d * d) };
				dg::vector::Vector<V> centreCircle(centre_ + h * (sphere.centre_ - centre_));
				dg::vector::Vector<V> normal(sphere.centre_ - centre_);
				dg::surf::Plane<V, V> planeCircle(normal, centreCircle);
				return planeCircle;
			}

			template<typename V, typename U,
				std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<U>, std::nullptr_t> = nullptr>
			double computeHeightIntersection_(const Plane<U, V>& plane)
			{
				dg::curve::Circumference<T, V, U> circle(intersectionCircle(plane).value());
				auto centreCircle(circle.getCentre());
				return radius_ - distance(centreCircle, centre_);
			}

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

			//getters 
			const dg::vector::Vector<T>& getCentre() const
			{
				return centre_;
			}
			
			T getRadius() const
			{
				return radius_;
			}

			T getDiameter() const
			{
				return 2 * radius_;
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
				dg::vector::Vector<T> pointOnSphere(centre_[0] + radius_ * sin(theta) * cos(phi),
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

			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr>
			double sphericalCapArea(V height) 
			{
				return 2.0 * dg::math::PI * radius_ * height;
			}

			template<typename V, typename U,
				std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<U>, std::nullptr_t> = nullptr>
			double sphericalCapArea(const Plane<U,V>& plane)
			{
				auto height{ computeHeightIntersection_(plane)};
				return sphericalCapArea(height);
				
			}

			//enclosed volume
			double volume() 
			{
				return 4.0 * dg::math::PI * radius_ * radius_ * radius_ / 3.0;
			}

			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr>
			double sphericalCapVolume(V height)
			{
				return dg::math::PI * height * height * (3.0 * radius_ - height) / 3.0;
			}

			template<typename V, typename U,
				std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<U>, std::nullptr_t> = nullptr>
			double sphericalCapVolume(const Plane<U, V>& plane) //will throw exception if no intersection
			{
				auto height{ computeHeightIntersection_(plane)};
				return sphericalCapVolume(height);

			}

			//antipodal point
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr>
			dg::vector::Vector<V> antipodal(const dg::vector::Vector<V>& point) //assumes point is on the sphere
			{
				dg::vector::Vector<V> antipodal(2.0 * centre_ - point);
				return antipodal;
			}
			//intersection with plane
			template<typename V, typename U,
				std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<U>, std::nullptr_t> = nullptr>
			bool hasIntersection(const Plane<V, U>& plane) const
			{
				auto signedDistance = plane.signedDistanceFrom(centre_);
				return fabs(signedDistance) < radius_;
			}

			template<typename V, typename U,
				std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<U>, std::nullptr_t> = nullptr>
			bool isTangent(const Plane<V, U>& plane) const 
			{
				return dg::math::almostEqualRelativeAndAbs(radius_, fabs(plane.signedDistanceFrom(centre_)), 1.0e-6);
			}

			template<typename V, typename U,
				std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<U>, std::nullptr_t> = nullptr>
			std::optional<dg::curve::Circumference<T,V,U>> intersectionCircle(const Plane<V, U>& plane) const
			{
				std::optional<dg::curve::Circumference<T,V,U>> intersectionCircle;
				if (hasIntersection(plane)) {
					auto signedDistance = plane.signedDistanceFrom(centre_);
					intersectionCircle = dg::curve::Circumference<T, V, U>(
						sqrt(radius_ * radius_ - signedDistance * signedDistance),
						Plane<V, U>(plane.getNormal(), dg::vector::Vector<U>(centre_ + signedDistance * plane.getNormal()))
						);
				}
				return intersectionCircle;
			}

			template<typename V, typename U,
				std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<U>, std::nullptr_t> = nullptr>
			std::optional<dg::vector::Vector<T>> tangent(const Plane<V, U>& plane) const
			{
				dg::vector::Vector<T> tangentPoint;
				if (isTangent(plane)) {
					dg::vector::Vector<V> normal(plane.getNormal());
					normal = ((angle(centre_, normal) < dg::math::PI) ? -1.0 : 1.0) * normal;
					tangentPoint = dg::vector::Vector<T>(centre_ + radius_ * normal);
				}
				return tangentPoint;
			}

			//intersection with other sphere
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr>
			bool hasIntersection(const Sphere<V>& sphere) const 
			{
				auto dist{ distance(centre_, sphere.centre_) };
				return (dist < radius_ + sphere.radius_) && (std::min(radius_, sphere.radius_) + dist > std::max(radius_, sphere.radius_));
			}

			template<typename V, 
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr>
			bool isTangent(const Sphere<V>& sphere) const
			{
				return isInternalTangent_<V>(sphere) || isExternalTangent_<V>(sphere);
			}

			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr>
			std::optional<dg::curve::Circumference<T, V, V>> intersectionCircle(const Sphere<V>& sphere) const
			{
				std::optional<dg::curve::Circumference<T, V, V>> intersectionCircle;
				if (hasIntersection(sphere))
				{
					auto d{ distance(centre_, sphere.centre_) };
					auto h{ 0.5 + (radius_ * radius_ - sphere.radius_ * sphere.radius_) / (2.0 * d * d) };
					auto radiusCircle{ sqrt(radius_ * radius_ - h * h * d * d) };
					intersectionCircle = dg::curve::Circumference<T, V, V>(radiusCircle, intersectionPlane_(sphere));
				}
				return intersectionCircle;
			}

			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr>
			std::optional<Plane<V,V>> intersectionPlane(const Sphere<V>& sphere) const
			{
				std::optional<Plane< V, V>> intersectionPlane;
				if (hasIntersection(sphere))
				{
					intersectionPlane = intersectionPlane_(sphere);
				}
				return intersectionPlane;

			}

			template<typename V, 
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr>
			std::optional<dg::vector::Vector<T>> tangent(const Sphere<V>& sphere) const
			{
				dg::vector::Vector<T> tangentPoint;
				if (isTangent(sphere)) {
					auto dist{ distance(centre_,sphere.centre_) };
					tangentPoint = dg::vector::Vector(centre_ + (radius_ / dist) * (sphere.centre_ - centre_));
				}
				return tangentPoint;
			}

			//great circle
			template<typename V, typename U,
				std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<U>, std::nullptr_t> = nullptr>
			dg::curve::Circumference<T, V, U> greatCircle(V theta, U phi) 
			{
				dg::vector::Vector<T> vectorAlongPlane(sin(theta) * cos(phi),
					sin(theta) * sin(phi),
					cos(phi));
				Plane<T, T> planeGreatCircle(vectorAlongPlane.perpendicular(), centre_);
				dg::curve::Circumference<T, V, U> geodesic(greatCircle(planeGreatCircle));
				return geodesic;
			}

			template<typename V, typename U,
				std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<U>, std::nullptr_t> = nullptr>
			dg::curve::Circumference<T, V, U> greatCircle(const Plane<U, V>& plane) 
			{
				dg::curve::Circumference<T, V, U> geodesic(radius_, plane);
				return geodesic; 
			}

			//Principal curvature
			auto principalCurvature() const
			{
				return -1.0 / radius_;
			}

			//Gaussian curvature
			auto guassianCurvature() const
			{
				return 1.0 / radius_ * radius_;
			}

			//Mean curvature
			auto meanCurvature() const
			{
				return -2.0 / radius_;
			}
		};

	} //namespace surf

} //namespace dg

#endif //!SPHERE_H


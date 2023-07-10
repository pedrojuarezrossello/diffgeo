#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#include <iostream>
#include <cmath>
#include <type_traits>
#include "../utils/type_traits.h"
#include "../utils/maths.h"
#include "../vector/vector_expression.h"


namespace dg {

	namespace vector {

		template<typename T>
		class Vector
		{
			static_assert(std::is_floating_point_v< T >, "Vector<T> only accepts numeric types!");
			using index = size_t;
			std::vector<T> vector_;

		public:

			Vector() = default;

			explicit Vector(std::vector<T>& v) : vector_(v) {}

			explicit Vector(std::initializer_list<T> initializer) : vector_(initializer) {}

			explicit Vector(T a, T b, T c) : vector_{ a,b,c } {}

			template<typename CallableObject, typename ... Args>
			Vector(const VectorExpression<CallableObject, Args...>& rhs) : vector_{ 0,0,0 }
			{
				for (size_t i = 0; i < 3; i++)
				{
					vector_[i] = rhs[i];
				}

			}

			template<typename CallableObject, typename ... Args>
			Vector& operator=(const VectorExpression<CallableObject, Args...>& rhs)
			{
				for (size_t i = 0; i < 3; i++)
				{
					this->vector_[i] = rhs[i];
				}

				return *this;
			}

			T operator[] (size_t i)
			{
				return vector_[i];
			}

			T operator[] (size_t i) const
			{
				return vector_[i];
			}

			friend std::ostream& operator<<(std::ostream& os, Vector const& vector)
			{
				os << "(";
				for (size_t i = 0; i < 2; ++i) {
					os << vector.vector_[i] << ",";
				}
				os << vector.vector_[2];
				os << ")";

				return os;
			}

			template<typename U,
				std::enable_if_t<std::is_floating_point_v<U>, std::nullptr_t> = nullptr >
			friend auto operator*(U lhs, Vector& rhs)
			{
				for (size_t i = 0; i < 3; i++) {
					rhs.vector_[i] *= lhs;
				}
				return rhs;
			}

			template<typename U,
				std::enable_if_t<std::is_floating_point_v<U>, std::nullptr_t> = nullptr >
			friend auto operator*(U lhs, const Vector& rhs)
			{
				Vector<T> scaledVector(rhs.vector_[0] * lhs, rhs.vector_[1] * lhs, rhs.vector_[2] * lhs);
				return scaledVector;
			}

			template<typename RHS>
			auto operator*(const RHS& rhs) const
			{
				static_assert(is_vector_or_expression_t<RHS>, "Can only dot by a vector or a vector expression");
				return vector_[0] * rhs[0] + vector_[1] * rhs[1] + vector_[2] * rhs[2];
			}

			template<typename U,
				std::enable_if_t<std::is_floating_point_v<U>, std::nullptr_t> = nullptr >
			friend double distance(const Vector& vec1, const Vector<U>& vec2)
			{
				const Vector<T> diffVec = vec1 - vec2;
				return diffVec.norm();
			}

			double norm() const
			{
				return sqrt(vector_[0] * vector_[0] + vector_[1] * vector_[1] + vector_[2] * vector_[2]);
			}

			template<typename Metric>
			auto norm(Metric metric_function) const
			{
				return metric_function(vector_[0], vector_[1], vector_[2]);
			}

			template<typename RHS>
			bool isPerpendicular(const RHS& vect)
			{
				return (dg::math::isZero(vect * (*this))) ? true : false;
			}

			template<typename RHS>
			Vector& orthogonalProjection(const RHS& vect)
			{
				if (dg::math::isZero(this->norm())) throw std::exception("Base vector can't be null");
				return (((*this) * vect) / pow(this->norm(), 2)) * (*this);
			}

			void normalise()
			{
				T norm = this->norm();

				if (dg::math::almostEqualRelativeAndAbs(norm, 0.0, 1.0e-6)) return;

				for (auto& component : vector_)
				{
					component *= 1.0 / norm;
				}
			}

			Vector perpendicular() const
			{
				Vector perpendicular(copysign(this->vector_[2], this->vector_[0]),
					copysign(this->vector_[2], this->vector_[1]),
					-copysign(abs(this->vector_[0]) + abs(this->vector_[1]), this->vector_[2]));

				perpendicular.normalise();

				return perpendicular; //(N)RVO
			}

			template<typename U,
				std::enable_if_t<std::is_floating_point_v<T>, std::nullptr_t> = nullptr>
			friend double angle(const Vector& vec1, const Vector<U>& vec2)
			{
				if (dg::math::isZero(vec1.norm()) || dg::math::isZero(vec2.norm())) throw std::exception("Vectors can't be null");
				return atan2(cross_product(vec1, vec2).norm(), vec1 * vec2)+dg::math::PI;
			}

			template<typename U,
				std::enable_if_t<std::is_floating_point_v<U>, std::nullptr_t> = nullptr>
			friend bool areParallel(const Vector& vec1, const Vector<U>& vec2)
			{
				return (dg::math::isZero(vec1.vector_[2] * vec2.vector_[0] - vec1.vector_[0] * vec2.vector_[2]) &&
					dg::math::isZero(vec1.vector_[1] * vec2.vector_[0] - vec1.vector_[0] * vec2.vector_[1]) &&
					dg::math::isZero(vec1.vector_[2] * vec2.vector_[1] - vec1.vector_[1] * vec2.vector_[2]));
			}

			bool operator==(const Vector& vec)
			{
				return dg::math::almostEqualRelativeAndAbs(this->vector_[0], vec.vector_[0], 1.0e-6) &&
					dg::math::almostEqualRelativeAndAbs(this->vector_[1], vec.vector_[1], 1.0e-6) &&
					dg::math::almostEqualRelativeAndAbs(this->vector_[2], vec.vector_[2], 1.0e-6);
			}
		};
	
	} //namespace vector

}//namespace dg

#endif //!VECTOR_H
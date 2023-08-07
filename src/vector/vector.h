#ifndef DG_VECTOR_H
#define DG_VECTOR_H
#include "../utils/type_traits.h"
#include "../utils/maths.h"
#include "../vector/vector_expression.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <type_traits>

namespace dg {

	namespace vector {

		template<typename T>
		class Vector
		{
			static_assert(std::is_floating_point_v< T >, "Vector<T> only accepts numeric types!");
			using index = size_t;
			std::vector<T> vector_;

			template<typename U>
			void scale_(U scalar)
			{
				for (auto& component : vector_) {
					component *= scalar;
				}
			}

			Vector perpendicular_()
			{
				return copysign(this->vector_[2], this->vector_[0]),
					copysign(this->vector_[2], this->vector_[1]),
					-copysign(abs(this->vector_[0]) + abs(this->vector_[1]), this->vector_[2]);
			}

		public:

			//Constructors

			Vector() = default;

			explicit Vector(std::vector<T>& v) : vector_(v) {}

			explicit Vector(std::initializer_list<T> initializer) : vector_(initializer) {}

			explicit Vector(T a, T b, T c) : vector_{ a,b,c } {}

			//Note these don't stop the compiler from generating copy/move constructors

			template<typename CallableObject, typename ... Args>
			Vector(const VectorExpression<CallableObject, Args...>& rhs) : vector_{ rhs[0],rhs[1],rhs[2]} {}

			template<typename CallableObject, typename ... Args>
			Vector& operator=(const VectorExpression<CallableObject, Args...>& rhs)
			{
				for (size_t i = 0; i < 3; i++)
				{
					this->vector_[i] = rhs[i];
				}

				return *this;
			}

			//Equality 

			bool operator==(const Vector& vec)
			{
				return dg::math::almostEqualRelativeAndAbs(this->vector_[0], vec.vector_[0], 1.0e-6) &&
					dg::math::almostEqualRelativeAndAbs(this->vector_[1], vec.vector_[1], 1.0e-6) &&
					dg::math::almostEqualRelativeAndAbs(this->vector_[2], vec.vector_[2], 1.0e-6);
			}

			//Component fetching

			T& operator[] (size_t i)
			{
				return vector_[i];
			}

			const T operator[] (size_t i) const
			{
				return vector_[i];
			}

			//Printing util

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

			//Scaling

			template<typename U,
				std::enable_if_t<std::is_floating_point_v<U>, std::nullptr_t> = nullptr >
			auto& operator*=(U rhs)
			{
				scale_(rhs);
				return *this;
			}

			template<typename U,
				std::enable_if_t<std::is_floating_point_v<U>, std::nullptr_t> = nullptr >
			friend auto operator*(U lhs, const Vector& rhs)
			{
				Vector<T> scaledVector(rhs.vector_[0] * lhs, rhs.vector_[1] * lhs, rhs.vector_[2] * lhs);
				return scaledVector;
			}

			//Dot product

			template<typename RHS>
			auto operator*(const RHS& rhs) const
			{
				return vector_[0] * rhs[0] + vector_[1] * rhs[1] + vector_[2] * rhs[2];
			}

			//Norms and metrics

			auto norm() const
			{
				return sqrt(vector_[0] * vector_[0] + vector_[1] * vector_[1] + vector_[2] * vector_[2]);
			}

			template<typename Metric>
			auto norm(Metric metric_function) const
			{
				return metric_function(vector_[0], vector_[1], vector_[2]);
			}

			void normalise()
			{
				T norm = this->norm();

				if (dg::math::almostEqualRelativeAndAbs(norm, 0.0, 1.0e-6)) return;

				scale_(1.0 / norm);
			}

			template<typename Metric>
			void normalise(Metric metric_function)
			{
				T norm = this->norm<Metric>(metric_function);

				if (dg::math::almostEqualRelativeAndAbs(norm, 0.0, 1.0e-6)) return;

				scale_(1.0 / norm);
			}


			//Perpendicular & parallel

			Vector perpendicular() const
			{
				auto perpVector(perpendicular_());
				perpVector.normalise();

				return perpVector; //(N)RVO
			}

			template<typename Metric>
			Vector perpendicular(Metric metric_function) const
			{
				auto perpVector(perpendicular_());
				perpVector.normalise<Metric>(metric_function);

				return perpVector; //(N)RVO
			}

			template<typename RHS>
			bool isPerpendicular(const RHS& vect)
			{
				return (dg::math::isZero(vect * (*this))) ? true : false;
			}

			template<typename U,
				std::enable_if_t<std::is_floating_point_v<U>, std::nullptr_t> = nullptr>
			friend bool areParallel(const Vector& vec1, const Vector<U>& vec2)
			{
				return (dg::math::isZero(vec1.vector_[2] * vec2.vector_[0] - vec1.vector_[0] * vec2.vector_[2]) &&
					dg::math::isZero(vec1.vector_[1] * vec2.vector_[0] - vec1.vector_[0] * vec2.vector_[1]) &&
					dg::math::isZero(vec1.vector_[2] * vec2.vector_[1] - vec1.vector_[1] * vec2.vector_[2]));
			}

			//Orthogonal projection

			template<typename RHS>
			Vector& orthogonalProjection(const RHS& vect)
			{
				return (((*this) * vect) / pow(this->norm(), 2)) * (*this);
			}

			//Angle

			template<typename U,
				std::enable_if_t<std::is_floating_point_v<T>, std::nullptr_t> = nullptr>
			friend double angle(const Vector& vec1, const Vector<U>& vec2)
			{
				return atan2(cross_product(vec1, vec2).norm(), vec1 * vec2)+dg::math::PI;
			}
		};

		//Distances
		template<typename U, typename LHS, typename RHS,
			std::enable_if_t<std::is_floating_point_v<U>, std::nullptr_t> = nullptr >
		auto distance(const LHS& vec1, const RHS& vec2)
		{
			const Vector<U> diffVec(vec1 - vec2);
			return diffVec.norm();
		}

		template<typename U, typename LHS, typename RHS, typename Metric,
			std::enable_if_t<std::is_floating_point_v<U>, std::nullptr_t> = nullptr >
		auto distance(const LHS& vec1, const RHS& vec2)
		{
			const Vector<U> diffVec(vec1 - vec2);
			return diffVec.norm<Metric>();
		}
	
	} //namespace vector

} //namespace dg

#endif //!DG_VECTOR_H
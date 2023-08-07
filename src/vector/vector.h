/** @file vector.h
 *	@ingroup Basic vectors
 *  @brief 3D vector class.
 * 
 *	It's the basis of the library, and it
 *	represents both a point in three dimensions 
 *	as well as a directed vector. It features the major 
 *	linear algebra related operations.
 * 
 *  @author Pedro Juarez Rossello (pedrojuarezrossello)
 *  @bug No known bugs.
 */

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

		/** 
		* @headerfile vector.h "src/vector/vector.h"
		* 
		* @tparam T Type of the vector components. 
		* It must be a floating-point type.
		*/
		template<typename T>
		class Vector
		{
			static_assert(std::is_floating_point_v< T >, "Vector<T> only accepts numeric types!");

			/**
			* @brief Underlying data structure is a three element array of type T.
			*/
			std::array<T,3> vector_;

			//Scaling util function
			template<typename U>
			void scale_(U scalar) 
			{
				for (auto& component : vector_) {
					component *= scalar;
				}
			}

			//Perpendicular util function
			Vector perpendicular_() const
			{
				return Vector(copysign(this->vector_[2], this->vector_[0]),
					copysign(this->vector_[2], this->vector_[1]),
					-copysign(abs(this->vector_[0]) + abs(this->vector_[1]), this->vector_[2]));
			}

		public:

			//Constructors

			/**
			* @brief Default constructor.
			*/
			Vector() = default;

			/**
			* @brief Constructor from `std::array`.
			* 
			* @param v 3 element array of type T.
			*/
			explicit Vector(std::array<T,3>& v) : vector_(v) {}

			/**
			* @brief Constructor from `std::initializer_list`.
			*
			* @param initializer 3 element initializer list of type T.
			*/
			explicit Vector(std::initializer_list<T> initializer) : vector_(initializer) {}

			/**
			* @brief Constructor from three ordered elements of type T.
			*
			* @param a First component of type T.
			* @param b Second component of type T.
			* @param c Third component of type T.
			*/
			explicit Vector(T a, T b, T c) : vector_{ a,b,c } {}

			//Note these don't stop the compiler from generating copy/move constructors

			/**
			* @brief Copy constructor from a VectorExpression.
			* 
			* @details It copies a given vector expression into a new vector,
			* for example v_1+2*v_2. Note that the expression isn't actually 
			* evaluated until it's forced to be turned into a vector.
			* 
			* @tparam CallableObject Operation to be perfomed on Args...
			* @tparam Args Argument pack of terms involved in the expression: v_1, v_2, etc.
			* @param rhs VectorExpression object (any linear combination of vectors).
			*/
			template<typename CallableObject, typename ... Args>
			Vector(const VectorExpression<CallableObject, Args...>& rhs) : vector_{ rhs[0],rhs[1],rhs[2]} {}

			/**
			* @brief Copy assignment constructor from a VectorExpression.
			*
			* @details It copies a given vector expression into an existing vector,
			* for example v_1+2*v_2. Note that the expression isn't actually
			* evaluated until it's forced to be turned into a vector.
			* 
			* @tparam CallableObject Operation to be perfomed on Args...
			* @tparam Args Argument pack of terms involved in the expression: v_1, v_2, etc.
			* @param rhs VectorExpression object (any linear combination of vectors).
			* 
			* @return Reference to the vector we're operating on.
			*/
			template<typename CallableObject, typename ... Args>
			Vector& operator=(const VectorExpression<CallableObject, Args...>& rhs)
			{
				for (size_t i = 0; i < 3; i++)
				{
					this->vector_[i] = rhs[i];
				}

				return *this;
			}

			//Equality operators

			/**
			* @build Equality operator between two vectors.
			* 
			* @details It compares two vectors componentwise by using
			* a relative difference measure for floating point types as opposed to
			* considering (naively) the absolute value of the difference.
			* 
			* @tparam U Type of the right hand side vector.
			* @param vec Vector of type U to be compared.
			* 
			* @return bool
			*/
			template<typename U,
				std::enable_if_t<std::is_floating_point_v<U>, std::nullptr_t> = nullptr >
			bool operator==(const Vector& vec)
			{
				return dg::math::almostEqualRelativeAndAbs(this->vector_[0], vec.vector_[0], 1.0e-11) &&
					dg::math::almostEqualRelativeAndAbs(this->vector_[1], vec.vector_[1], 1.0e-11) &&
					dg::math::almostEqualRelativeAndAbs(this->vector_[2], vec.vector_[2], 1.0e-11);
			}

			/**
			* @build Equality operator between a vector and a vector expression.
			*
			* @details It compares a vector and a vector expression componentwise 
			* (and thus forcing the evaluation of the vector expression) by using
			* a relative difference measure for floating point types as opposed to
			* considering (naively) the absolute value of the difference.
			*
			* @tparam CallableObject Operation to be perfomed on Args...
			* @tparam Args Argument pack of terms involved in the expression: v_1, v_2, etc.
			* @param rhs VectorExpression object (any linear combination of vectors).
			* 
			* @return bool
			*/
			template<typename CallableObject, typename ... Args>
			bool operator==(const VectorExpression<CallableObject, Args...>& rhs) 
			{
				return dg::math::almostEqualRelativeAndAbs(this->vector_[0], rhs.vector_[0], 1.0e-11) &&
					dg::math::almostEqualRelativeAndAbs(this->vector_[1], rhs.vector_[1], 1.0e-11) &&
					dg::math::almostEqualRelativeAndAbs(this->vector_[2], rhs.vector_[2], 1.0e-11);
			}

			//Component fetching

			/** 
			* @brief Component fetching (non-const)
			* 
			* @detail It fetches the ith component and returns a reference
			* to this element, allowing us to modify the component - read and write.
			* 
			* @param i Index we'd like to fetch.
			* 
			* @return A reference to the component.
			*/
			T& operator[] (size_t i)
			{
				return vector_[i];
			}

			/**
			* @brief Component fetching (const)
			*
			* @detail It fetches the ith component and returns 
			* a const copy of this element, so it's read only.
			*
			* @param i Index we'd like to fetch.
			*
			* @return A const copy the component.
			*/
			const T operator[] (size_t i) const
			{
				//not const T& as T is a built-in type
				return vector_[i]; 
			}

			//Printing util

			/**
			* @brief Printing util for Vector class
			* 
			* @details It prints a vector in the format \f$(x,y,z)\f$.
			*/
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

			/**
			* @brief Scaling (in-place) of a vector
			* 
			* @details It multiplies a vector componentwise
			* by a given scalar rhs. It modifies the vector.
			* 
			* @tparam U Type of multiplying scalar.
			* @param rhs Scalar of type U.
			* 
			* @return The original vector, modified.
			*/
			template<typename U,
				std::enable_if_t<std::is_floating_point_v<U>, std::nullptr_t> = nullptr >
			auto& operator*=(U rhs)
			{
				scale_(rhs);
				return *this;
			}

			/**
			* @brief Scaling - new vector.
			*
			* @details It multiplies a vector componentwise
			* by a given scalar rhs and returns the result in a newly created Vector.
			*
			* @tparam U Type of multiplying scalar.
			* @param rhs Scalar of type U.
			*
			* @return New vector that's equal to the scaled original vector.
			*/
			template<typename U,
				std::enable_if_t<std::is_floating_point_v<U>, std::nullptr_t> = nullptr >
			friend auto operator*(U lhs, const Vector& rhs)
			{
				Vector<T> scaledVector(rhs.vector_[0] * lhs, rhs.vector_[1] * lhs, rhs.vector_[2] * lhs);
				return scaledVector;
			}

			//Dot product

			/**
			* @brief Dot product.
			* 
			* @details It works out the dot product of two Vector objects,
			* or of a Vector object and a VectorExpression object.
			* 
			* @tparam RHS It can be a Vector or a VectorExpression.
			* @param rhs Vector or VectorExpression object 
			* to be dotted against the working vector.
			* 
			* @return A scalar equal to the dot product.
			*/
			template<typename RHS>
			auto operator*(const RHS& rhs) const
			{
				return vector_[0] * rhs[0] + vector_[1] * rhs[1] + vector_[2] * rhs[2];
			}

			//Norms and metrics

			/**
			* @brief Euclidean norm.
			*
			* @details It works out Euclidean norm i.e. the length of a vector.
			* 
			* @return A scalar equal to the norm.
			*/
			auto norm() const
			{
				return sqrt(vector_[0] * vector_[0] + vector_[1] * vector_[1] + vector_[2] * vector_[2]);
			}

			/**
			* @brief General norm.
			*
			* @details It works out the norm under a different given metric.
			*
			* @tparam Metric Any function type having signature F (F,F,F), 
			* where F is any floating point type.
			* @param metric_function A function representing the metric.
			* 
			* @return A scalar equal to the norm.
			*/
			template<typename Metric>
			auto norm(Metric metric_function) const
			{
				return metric_function(vector_[0], vector_[1], vector_[2]);
			}

			/**
			* @brief Normalise a vector.
			* 
			* @details It modifies a vector to have length 1.
			* 
			* @return void.
			*/
			void normalise()
			{
				T norm = this->norm();

				if (dg::math::almostEqualRelativeAndAbs(norm, 0.0, 1.0e-6)) return;

				scale_(1.0 / norm);
			}

			/**
			* @brief Normalise a vector.
			*
			* @details It modifies a vector to have length 1 
			* under a given metric.
			*
			* @tparam Metric Any function type having signature F (F,F,F), 
			* where F is any floating point type.
			* @param metric_function A function representing the metric.
			* 
			* @return void.
			*/
			template<typename Metric>
			void normalise(Metric metric_function)
			{
				T norm = this->norm<Metric>(metric_function);

				if (dg::math::almostEqualRelativeAndAbs(norm, 0.0, 1.0e-6)) return;

				scale_(1.0 / norm);
			}


			//Perpendicular & parallel

			/**
			* @brief It finds a unit perpendicular.
			* 
			* @details Note the perpendicular is not unique.
			* 
			* @return Perpendicular Vector object.
			*/
			Vector perpendicular() const
			{
				auto perpVector(perpendicular_());
				perpVector.normalise();

				return perpVector; //(N)RVO
			}

			/**
			* @brief It finds a unit perpendicular in a given metric.
			*
			* @details Note the perpendicular is not unique.
			*
			* @tparam Metric Any function type having signature F (F,F,F), 
			* where F is any floating point type.
			* @param metric_function A function representing the metric.
			* 
			* @return Perpendicular Vector object.
			*/
			template<typename Metric>
			Vector perpendicular(Metric metric_function) const
			{
				auto perpVector(perpendicular_());
				perpVector.normalise<Metric>(metric_function);

				return perpVector; //(N)RVO
			}

			/**
			* @brief Test of perpendicularity.
			* 
			* @details It finds out if a given Vector or VectorExpression is 
			* perpendicular to the working Vector instance.
			* 
			* @tparam RHS It can be a Vector or a VectorExpression.
			* @param rhs Vector or VectorExpression object 
			* to be compared against the working vector.
			* 
			* @return bool
			*/
			template<typename RHS>
			bool isPerpendicular(const RHS& vect)
			{
				return (dg::math::isZero(vect * (*this))) ? true : false;
			}

			/**
			* @brief Test of parallel vectors.
			*
			* @details It finds out if a given Vector is
			* parallel to the working Vector instance by computing their cross product.
			*
			* @tparam U Type of the right hand side vector.
			* @param vec Vector of type U to be compared.
			* 
			* @return bool
			*/
			template<typename U,
				std::enable_if_t<std::is_floating_point_v<U>, std::nullptr_t> = nullptr>
			friend bool areParallel(const Vector& vec1, const Vector<U>& vec2)
			{
				return (dg::math::isZero(vec1.vector_[2] * vec2.vector_[0] - vec1.vector_[0] * vec2.vector_[2]) &&
					dg::math::isZero(vec1.vector_[1] * vec2.vector_[0] - vec1.vector_[0] * vec2.vector_[1]) &&
					dg::math::isZero(vec1.vector_[2] * vec2.vector_[1] - vec1.vector_[1] * vec2.vector_[2]));
			}

			//Orthogonal projection

			/**
			* @brief Compute orthogonal projections.
			* 
			* @details It computes the orthogonal projection of the given Vector or 
			* VectorExpression on the working vector (the on referred by the this pointer).
			* 
			* @tparam RHS It can be a Vector or a VectorExpression.
			* @param rhs Vector or VectorExpression object 
			* to be compared against the working vector.
			* 
			* @return Vector object containing the orthogonal projection.
			*/
			template<typename RHS>
			Vector orthogonalProjection(const RHS& vect)
			{
				Vector proj((((*this) * vect) / pow(this->norm(), 2)) * (*this));
				return proj;
			}

			//Angle

			/**
			* @brief It works out the angle between two vectors.
			* 
			* @details The angle is between 0 and \f$ 2\pi \f$.
			* 
			* @tparam U Type of the second Vector (in argument list).
			* @param vec1 First vector.
			* @param vec2 Second vector.
			* 
			* @return Scalar.
			*/
			template<typename U,
				std::enable_if_t<std::is_floating_point_v<T>, std::nullptr_t> = nullptr>
			friend double angle(const Vector& vec1, const Vector<U>& vec2)
			{
				return atan2(cross_product(vec1, vec2).norm(), vec1 * vec2)+dg::math::PI;
			}
		};

		//Distances

			/**
			* @brief Euclidean distance.
			*
			* @details It works out Euclidean distnce between two
			* Vector or VectorExpression objects by computing the 
			* Euclidean norm of its difference.
			* 
			* @tparam U Floating point type we'd like the solution in.
			* @tparam LHS Operator 1 type.
			* @tparam RHS Operator 2 type.
			* @param vec1 Operator 1.
			* @param vec2 Operator 2.
			* 
			* @return A scalar equal to the distance.
			*/
		template<typename U, typename LHS, typename RHS,
			std::enable_if_t<std::is_floating_point_v<U>, std::nullptr_t> = nullptr >
		auto distance(const LHS& vec1, const RHS& vec2)
		{
			const Vector<U> diffVec(vec1 - vec2);
			return diffVec.norm();
		}

			/**
			* @brief General distance.
			*
			* @details It works out distnce in a given metric between two
			* Vector or VectorExpression objects by computing the norm of its difference.
			*
			* @tparam U Floating point type we'd like the solution in.
			* @tparam LHS Operator 1 type.
			* @tparam RHS Operator 2 type.
			* @tparam Metric Any function type having signature F (F,F,F), 
			* where F is any floating point type.
			* @param vec1 Operator 1.
			* @param vec2 Operator 2.
			* @param metric_function A function representing the metric.
			* 
			* @return A scalar equal to the distance.
			*/
		template<typename U, typename LHS, typename RHS, typename Metric,
			std::enable_if_t<std::is_floating_point_v<U>, std::nullptr_t> = nullptr >
		auto distance(const LHS& vec1, const RHS& vec2, Metric metric_function)
		{
			const Vector<U> diffVec(vec1 - vec2);
			return diffVec.norm<Metric>(metric_function);
		}
	
	} //namespace vector

} //namespace dg

#endif //!DG_VECTOR_H
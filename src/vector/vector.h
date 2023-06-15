#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#include <iostream>
#include <cmath>
#include <type_traits>
#include "../utils/type_traits.h"
#include "../utils/maths.h"
using std::vector;
using std::cout;
using std::ostream;


template<typename T>
class Vector
{
	static_assert(std::is_floating_point_v< T >, "Vector<T> only accepts numeric types!");
	using index = size_t;
	vector<T> vector_;

public:

	Vector() : vector_{0,0,0} {}

	explicit Vector(vector<T>& v) : vector_(v) {}

	explicit Vector(std::initializer_list<T> initializer) : vector_(initializer) {}

	explicit Vector(T a, T b, T c) : vector_{ a,b,c } {}

	template<typename CallableObject, typename ... Args>
	Vector(const VectorExpression<CallableObject, Args...>& rhs) : vector_{0,0,0}
	{
		for (size_t i =0 ; i<3;i++)
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

	friend ostream& operator<<(ostream& os, Vector const& vector)
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
	 friend auto&& operator*(U lhs, Vector& rhs)
	 {
		 for (size_t i = 0; i < 3; i++) {
			 rhs.vector_[i] *= lhs;
		 }
		 return rhs;
	 }

	template<typename RHS>
	auto operator*(const RHS& rhs) const;

	template<typename RHS>
	auto operator*(RHS&& rhs) const;

	template<typename U,
		std::enable_if_t<std::is_floating_point_v<U>, std::nullptr_t> = nullptr >
	friend double distance(const Vector& vec1, const Vector<U>& vec2)
	{
		const Vector<T> diffVec = vec1 - vec2;
		return diffVec.norm();
	}

	double norm() const;

	template<typename Metric>
	auto norm(Metric metric_function) const;

	template<typename RHS>
	bool isPerpendicular(const RHS& vect)
	{
		return (isZero(vect * (*this))) ? true : false;
	}

	template<typename RHS>
	Vector& orthogonalProjection(const RHS& vect)
	{
		if (isZero(this->norm())) throw std::exception("Base vector can't be null");
		return (((*this) * vect) / pow(this->norm(),2)) * (*this);
	}

	void normalise()
	{
		T norm = this->norm();

		if (almostEqualRelativeAndAbs(norm, 0.0, 1.0e-6)) return;

		for (auto& component : vector_)
		{
			component *= 1.0/norm;
		}
	}

	Vector perpendicular()
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
		if (isZero(vec1.norm()) || isZero(vec2.norm())) throw std::exception("Vectors can't be null");
		return atan2(cross_product(vec1, vec2).norm(), vec1*vec2);
	}

	template<typename U,
		std::enable_if_t<std::is_floating_point_v<T>, std::nullptr_t> = nullptr>
	friend bool areParallel(const Vector& vec1, const Vector<U>& vec2)
	{
		return (isZero(vec1.vector_[2] * vec2.vector_[0] - vec1.vector_[0] * vec2.vector_[2]) &&
			isZero(vec1.vector_[1] * vec2.vector_[0] - vec1.vector_[0] * vec2.vector_[1]) &&
			isZero(vec1.vector_[2] * vec2.vector_[1] - vec1.vector_[1] * vec2.vector_[2]));
	}

	template<typename U,
		std::enable_if_t<std::is_floating_point_v<T>, std::nullptr_t> = nullptr >
	bool equals(const Vector<U>& vec)
	{
		return almostEqualRelativeAndAbs(this->vector_[0], vec.vector_[0], 1.0e-6) &&
			almostEqualRelativeAndAbs(this->vector_[1], vec.vector_[1], 1.0e-6) &&
			almostEqualRelativeAndAbs(this->vector_[2], vec.vector_[2], 1.0e-6);
	}
};

template<typename T>
double Vector<T>::norm() const
{
	return sqrt(pow(static_cast<double>(vector_[0]), 2) + pow(static_cast<double>(vector_[1]), 2) + pow(static_cast<double>(vector_[2]), 2));
}

template<typename T>
template<typename Metric>
auto Vector<T>::norm(Metric metric_function) const
{
	return metric_function(vector_[0], vector_[1], vector_[2]);
}

template <typename T>
template <typename RHS>
auto Vector<T>::operator*(const RHS& rhs) const
{
	static_assert(is_vector_or_expression_t<RHS>, "Can only dot by a vector or a vector expression");
	return vector_[0] * rhs[0] + vector_[1] * rhs[1] + vector_[2] * rhs[2];
}

template <typename T>
template <typename RHS>
auto Vector<T>::operator*(RHS&& rhs) const
{
	static_assert(is_vector_or_expression_t<RHS>, "Can only dot by a vector or a vector expression");
	return vector_[0] * rhs[0] + vector_[1] * rhs[1] + vector_[2] * rhs[2];
}


#endif //!VECTOR_H
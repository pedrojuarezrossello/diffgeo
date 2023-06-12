#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#include <iostream>
#include <cmath>
#include <type_traits>
#include "../utils/type_traits.h"

using std::vector;
using std::cout;
using std::ostream;


template<typename T>
class Vector
{
	static_assert(std::is_arithmetic_v < T >, "Vector<T> only accepts numeric types!");
	using index = size_t;
	vector<T> vector_;

public:

	explicit Vector(vector<T>& v) : vector_(v) {}

	explicit Vector(std::initializer_list<T> init) : vector_(init) {}

	explicit Vector(T a, T b, T c) : vector_{ a,b,c } {}

	template <typename RHS>
	Vector& operator=(RHS const& rhs);

	auto operator[] (size_t i)
	{
		return vector_[i];
	}

	auto operator[] (size_t i) const
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
		os << ")\n";

		return os;
	}

	template<typename RHS>
	auto operator*(const RHS& rhs) const;

	double norm() const;

	template<typename Metric>
	auto norm(Metric metric_function) const;
};


template <typename T>
template <typename RHS>
Vector<T>& Vector<T>::operator=(RHS const& rhs)
{
	static_assert(is_vector_or_expression_t<RHS>, "Can only assign a vector or a vector expression");
	for (index i=0; i<3; i++)
	{
		vector_[i] = rhs[i];
	}

	return *this;
}

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


#endif //!VECTOR_H
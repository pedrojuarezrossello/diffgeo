#ifndef VECTOR_H
#define VECTOR_H
#include <vector>
#include <iostream>

using std::vector;
using std::cout;
using std::ostream;

template<typename T>
class Vector
{
	using index = size_t;
	vector<T> vector_;

public:

	explicit Vector(vector<T>& v) : vector_(v) {}

	explicit Vector(std::initializer_list<T> init) : vector_(init) {}

	explicit Vector(T a, T b, T c) : vector_{ a,b,c } {}

	template <typename RHS>
	Vector& operator=(RHS const& rhs);

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

	auto operator[] (size_t i)
	{
		return vector_[i];
	}

	auto operator[] (size_t i) const
	{

		return vector_[i];
	}
};

template <typename T>
template <typename RHS>
Vector<T>& Vector<T>::operator=(RHS const& rhs)
{
	for (index i=0; i<3; i++)
	{
		vector_[i] = rhs[i];
	}

	return *this;
}


#endif //!VECTOR_H
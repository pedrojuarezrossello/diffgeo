#ifndef OPERATIONS_H
#define OPERATIONS_H
#include "vector_expression.h"

template<typename LHS, typename RHS>
auto operator+(const LHS& lhs, const RHS& rhs)
{
	auto expression = VectorExpression{				//class template argument deduction C++17
		[](auto const& left, auto const& right)
		{
			return left + right;
		},lhs,rhs };

	return expression; //(N)RVO
}

template<typename LHS, typename RHS>
auto operator-(const LHS& lhs, const RHS& rhs)
{
	auto expression = VectorExpression{				//class template argument deduction C++17
		[](auto const& left, auto const& right)
		{
			return left - right;
		},lhs,rhs };

	return expression; //(N)RVO
}

//only for matrices and stuff
/*template<typename LHS, typename RHS>
auto operator*(const LHS& lhs, const RHS& rhs)
{
	auto expression = VectorExpression{				//class template argument deduction C++17
		[](auto const& left, auto const& right)
		{
			return left * right;
		},lhs,rhs };

	return expression; //(N)RVO
}*/




#endif //!OPERATIONS_H

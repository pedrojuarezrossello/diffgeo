#ifndef OPERATIONS_H
#define OPERATIONS_H
#include <type_traits>
#include "vector_expression.h"
#include "../utils/type_traits.h"

namespace dg {

	namespace vector {

		template<typename LHS, typename RHS,
			std::enable_if_t<is_vector_or_expression_t<RHS>&& is_vector_or_expression_t<LHS>, std::nullptr_t> = nullptr>
		auto operator+(const LHS& lhs, const RHS& rhs)
		{
			auto expression = VectorExpression{				//class template argument deduction C++17
				[](auto const& left, auto const& right)
				{
					return left + right;
				},lhs,rhs };

			return expression; //(N)RVO
		}

		template<typename LHS, typename RHS,
			std::enable_if_t<is_vector_or_expression_t<RHS>&& is_vector_or_expression_t<LHS>, std::nullptr_t> = nullptr>
		auto operator-(const LHS& lhs, const RHS& rhs)
		{
			auto expression = VectorExpression{				//class template argument deduction C++17
				[](auto const& left, auto const& right)
				{
					return left - right;
				},lhs,rhs };

			return expression; //(N)RVO
		}

		template<typename T, typename CallableObject, typename... Args,
			std::enable_if_t<std::is_floating_point_v<T>, std::nullptr_t> = nullptr>
		auto operator*(T lhs, const VectorExpression<CallableObject, Args...>& rhs)
		{
			Vector<T> vec(rhs[0], rhs[1], rhs[2]);

			return lhs * vec;
		}

		template<typename T, typename CallableObject, typename... Args,
			std::enable_if_t<std::is_floating_point_v<T>, std::nullptr_t> = nullptr>
		auto operator*(T lhs, VectorExpression<CallableObject, Args...>&& rhs)
		{
			Vector<T> vec(rhs[0], rhs[1], rhs[2]);

			return lhs * vec;
		}

		template<typename LHS, typename RHS,
			std::enable_if_t<is_vector_or_expression_t<RHS>&& is_vector_or_expression_t<LHS>, std::nullptr_t> = nullptr>
		auto cross_product(const LHS& lhs, const RHS& rhs)
		{
			auto cross = Vector(lhs[1] * rhs[2] - lhs[2] * rhs[1],
				lhs[2] * rhs[0] - lhs[0] * rhs[2],
				lhs[0] * rhs[1] - lhs[1] * rhs[0]);

			return cross; //(N)RVO
		}

		template<typename T, typename LHS, typename RHS,
			std::enable_if_t<is_vector_or_expression_t<RHS>&& is_vector_or_expression_t<LHS> && std::is_floating_point_v<T>, std::nullptr_t> = nullptr>
		Vector<T> midPoint(const LHS& lhs, const RHS& rhs)
		{
			Vector<T> mid(0.5 * (lhs + rhs));
			return mid;
		}

	} //namespace vector

} //namespace dg

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

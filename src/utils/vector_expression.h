#ifndef VECTOR_EXPRESSION_H
#define VECTOR_EXPRESSION_H
#include <tuple>
#include <iostream>
#include "../utils/type_traits.h"

using std::tuple;
using std::apply;

template <typename CallableObject, typename... Args>
class VectorExpression
{
	static_assert(std::conjunction_v<is_vector_or_expression<Args>...>,
		"Can only build a vector expression out of vectors and vector expressions");

	using index = size_t;
	tuple<Args const&...> operands_;
	CallableObject func_;

public:
	explicit VectorExpression(CallableObject func, Args const&... operands)
		: operands_(operands...), func_(func) {}

	auto operator[] (index const i) const
	{
		auto call_at =
			[this, &i](Args const&... operands)
		{
			return func_(operands[i]...);
		};

		return apply(call_at, operands_);
	}

	template<typename RHS>
	auto operator*(const RHS& rhs) const;
};

template<typename CallableObject, typename... Args>
template<typename RHS>
auto VectorExpression<CallableObject,Args...>::operator*(const RHS& rhs) const
{
	static_assert(is_vector_or_expression_t<RHS>, "Can only dot by a vector or a vector expression");
	return (*this)[0] * rhs[0] + (*this)[1] * rhs[1] + (*this)[2] * rhs[2];
}

#endif //!VECTOR_EXPRESSION_h


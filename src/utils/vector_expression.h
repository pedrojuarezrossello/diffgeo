#ifndef VECTOR_EXPRESSION_H
#define VECTOR_EXPRESSION_H
#include <tuple>
#include <iostream>

using std::tuple;
using std::apply;

template <typename CallableObject, class... Args>
class VectorExpression
{
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


};


#endif //!VECTOR_EXPRESSION_h


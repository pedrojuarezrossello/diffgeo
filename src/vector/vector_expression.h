#ifndef DG_VECTOR_EXPRESSION_H
#define DG_VECTOR_EXPRESSION_H
#include "../utils/type_traits.h"
#include <tuple>
#include <iostream>

namespace dg {

	namespace vector {

		template <typename CallableObject, typename... Args>
		class VectorExpression
		{
			static_assert(std::conjunction_v<is_vector_or_expression<Args>...>,
				"Can only build a vector expression out of vectors and vector expressions");

			using index = size_t;
			std::tuple<Args const&...> operands_;
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

				return std::apply(call_at, operands_);
			}

			template<typename RHS>
			auto operator*(const RHS& rhs) const;

			friend std::ostream& operator<<(std::ostream& os, VectorExpression const& vector)
			{
				os << "(";
				for (size_t i = 0; i < 2; ++i) {
					os << vector[i] << ",";
				}
				os << vector[2];
				os << ")\n";

				return os;
			}
		};

		template<typename CallableObject, typename... Args>
		template<typename RHS>
		auto VectorExpression<CallableObject, Args...>::operator*(const RHS& rhs) const
		{
			static_assert(is_vector_or_expression_t<RHS>, "Can only dot a vector or a vector expression");
			return (*this)[0] * rhs[0] + (*this)[1] * rhs[1] + (*this)[2] * rhs[2];
		}

	} //namespace vector

} //namespace dg

#endif //!DG_VECTOR_EXPRESSION_H


#ifndef REAL_FUNCTION_H
#define REAL_FUNCTION_H
#include <type_traits>
#include <functional>
#include <boost/math/differentiation/autodiff.hpp>

namespace dg {

	namespace math {

		template<typename T>
		class Function
		{
			using Funct = auto (T) -> T;

			std::function<T(T)> underlying_function;

			template<int Order, typename U,
				std::enable_if_t<std::is_floating_point_v<U>, std::nullptr_t> = nullptr>
			auto find_derivatives_helper_(U value)
			{
				auto const epsilonedValue = boost::math::differentiation::make_fvar<U, Order>(value);
				auto derivatives = underlying_function(epsilonedValue);
				return derivatives;
			}

		public:

			Function(std::function<Funct> func) : underlying_function(func) {}

			template<typename U,
				std::enable_if_t<std::is_floating_point_v<U>, std::nullptr_t> = nullptr>
			U operator()(U var)
			{
				auto derivatives = find_derivatives_helper_<0>(var);
				return derivatives.derivative(0);
			}

			template<int Order, typename U,
				std::enable_if_t<std::is_floating_point_v<U>, std::nullptr_t> = nullptr>
			U derivative(U var)
			{
				auto derivatives = find_derivatives_helper_<Order>(var);
				return derivatives.derivative(Order);
			}


		};

	} //namespace math

} //namespace dg



#endif //!REAL_FUNCTION_H


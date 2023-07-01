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

			std::function<Funct> underlying_function;

			template<int Order, typename U,
				std::enable_if_t<std::is_floating_point_v<U>, std::nullptr_t> = nullptr >
			auto find_derivatives_helper_(U value) const
			{
				auto const epsilonedValue = boost::math::differentiation::make_fvar<U, Order>(value);
				auto derivatives = underlying_function(epsilonedValue);
				return derivatives;
			}

		public:

			Function() {}

			Function(std::function<Funct> func) : underlying_function(func) {}

			template<typename U>
			U operator()(U var) 
			{
				return underlying_function(var).derivative(0);
			}

			T operator()(T var) const
			{
				return underlying_function(var);
			}

			template<int Order, typename U,
				std::enable_if_t<std::is_floating_point_v<U>, std::nullptr_t> = nullptr >
			U derivative(U var) const
			{
				auto derivatives = find_derivatives_helper_<Order>(var);
				return derivatives.derivative(Order);
			}

			template<int Order>
			auto derivative(T var) const
			{
				auto derivatives = underlying_function(var);
				return derivatives.derivative(Order);
			}

		};


	} //namespace math

} //namespace dg



#endif //!REAL_FUNCTION_H


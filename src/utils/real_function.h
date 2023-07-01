#ifndef REAL_FUNCTION_H
#define REAL_FUNCTION_H
#include <type_traits>
#include <functional>
#include <boost/math/differentiation/autodiff.hpp>
#include <iostream>

namespace dg {

	namespace math {

		template<typename _Ret, typename... Components>
		class Function
		{
			using Funct = auto (Components...) -> _Ret;

			std::function<Funct> underlying_function;
		
			template<int Order1, typename U>
			auto find_derivatives_helper_(U value) const
			{
				auto const epsilonedValue = boost::math::differentiation::make_fvar<U, Order1>(value);
				auto derivatives = underlying_function(epsilonedValue);
				std::cout << " I'm single valued";
				return derivatives;
			}

			

		public:

			template<int Order1, int Order2, typename ... Us>
			auto find_derivatives_helper_(Us ... values) const
			{
				using Tuple = std::tuple<Us...>;
				auto const epsilonedValues = boost::math::differentiation::make_ftuple<std::tuple_element<0,Tuple>::type, Order1, Order2>(values...);
				auto const& x = std::get<0>(epsilonedValues);
				auto const& y = std::get<1>(epsilonedValues);
				auto const derivatives = underlying_function(x,y);
				//std::cout << x << y;
				return derivatives;
			}

			Function() {}

			Function(std::function<Funct> func) : underlying_function(func) {}

			auto operator()(Components... var) const
			{
				return underlying_function(var...);
			}

			/*template<typename... U>
			auto operator()(U... var) 
			{
				return underlying_function(var...).derivative(0);
			}

			auto operator()(T... var) const
			{
				return underlying_function(var...);
			}

			template<int Order, typename... U,
				std::enable_if_t<std::is_floating_point_v<U...>, std::nullptr_t> = nullptr >
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
			*/
		};


	} //namespace math

} //namespace dg



#endif //!REAL_FUNCTION_H


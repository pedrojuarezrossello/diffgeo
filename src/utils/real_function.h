/**
* @file real_function.h
* @ingroup Maths
* @brief Real function in one or more variables.
* 
* @author Pedro Juarez Rossello (pedrojuarezrossello)
* @bug No known bugs.
*/

#ifndef DG_REAL_FUNCTION_H
#define DG_REAL_FUNCTION_H
#include <boost/math/differentiation/autodiff.hpp>
#include <type_traits>
#include <functional>

namespace dg {

	namespace math {

		/**
		* @headerfile real_function.h "src/utils/real_function.h"
		* 
		* @details The Function class plays a pivotal role
		* in the library as it allows us to very accurately differentiate
		* functions in one or more variables using the autodifferentiation 
		* library in Boost.Math. The price to pay is a slightly awkward syntax.
		* 
		*/
		template<typename _Ret, typename... Components>
		class Function
		{
			using Funct = auto (Components...) -> _Ret;

			//The underlying function is std::function
			std::function<Funct> underlying_function;
		
			//Helper method to single-variable derivatives
			template<int Order1, typename U>
			auto find_derivatives_helper_(U value) const
			{
				auto const epsilonedValue = boost::math::differentiation::make_fvar<U, Order1>(value);
				auto derivatives = underlying_function(epsilonedValue);
				return derivatives;
			}

			//Helper method to find partial derivatives
			template<int Order1, int Order2, typename ... Us>
			auto find_derivatives_helper_(Us ... values) const
			{
				using Tuple = std::tuple<Us...>;
				auto const epsilonedValues = boost::math::differentiation::make_ftuple<std::tuple_element<0, Tuple>::type, Order1, Order2>(values...);
				auto const& x = std::get<0>(epsilonedValues);
				auto const& y = std::get<1>(epsilonedValues);
				auto const derivatives = underlying_function(x, y);
				return derivatives;
			}

			//Derivative function util
			template<int Order, typename U,
				std::enable_if_t<std::is_floating_point_v<U>, std::nullptr_t> = nullptr >
			U derivative_(U var) const
			{
				auto derivatives = find_derivatives_helper_<Order>(var);
				return derivatives.derivative(Order);
			}

		public:

			/**
			* @brief Default constructor
			*/
			Function() = default;

			Function(std::function<Funct> func) : underlying_function(func) {}

			auto operator()(Components... var) const
			{
				return underlying_function(var...);
			}

			template<typename... U>
			auto operator()(U... var) 
			{
				return underlying_function(var...).derivative(static_cast<int>(var-var)...);
			}

			template<int Order, typename U>
			U derivative(U var) const
			{
				return derivative_<Order, U>(var);
			}

			template<int Order, int Order2, typename... Us>
			auto derivative(Us... var) const
			{
				auto derivatives = find_derivatives_helper_<Order,Order2>(var...);
				return derivatives.derivative(Order,Order2);
			}

			template<int Order>
			auto derivative(Components... var) const
			{
				auto derivatives = underlying_function(var...);
				return derivatives.derivative(Order);
			}
			
			template<int Order1, int Order2>
			auto derivative(Components... var) const
			{
				auto derivatives = underlying_function(var...);
				return derivatives.derivative(Order1,Order2);
			}
		};


	} //namespace math

} //namespace dg



#endif //!DG_REAL_FUNCTION_H


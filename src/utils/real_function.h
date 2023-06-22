#ifndef REAL_FUNCTION_H
#define REAL_FUNCTION_H
#include <type_traits>
#include <functional>

namespace dg {

	namespace math {

		template<typename T,
			std::enable_if_t<std::is_floating_point_v<T>, std::nullptr_t> = nullptr>
		struct Function
		{
			using Funct = auto (double) -> double;

			Function(std::function<Funct> func) : underlying_function(func) {}

			std::function<T(T)> underlying_function;
			T operator()(T var)
			{
				return underlying_function(var);
			}

		};

	} //namespace math

} //namespace dg



#endif //!REAL_FUNCTION_H


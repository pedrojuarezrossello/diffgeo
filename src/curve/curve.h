#ifndef CURVE_H
#define CURVE_H
#include <type_traits>
#include "../utils/real_function.h"

namespace dg {

	namespace curve {

		template<typename T,
			std::enable_if_t<std::is_floating_point_v<T>, std::nullptr_t> = nullptr>
		class Curve
		{
			dg::math::Function<T> X_;
			dg::math::Function<T> Y_;
			dg::math::Function<T> Z_;

		public:
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			dg::vector::Vector<T> operator()(V var)
			{
				dg::vector::Vector<T> evaluation(X_(var), Y_(var), Z_(var));
				return evaluation;
			}

			Curve() = default;

			Curve(const dg::math::Function<T>& X, const dg::math::Function<T>& Y, const dg::math::Function<T>& Z)
				: X_(X), Y_(Y), Z_(Z)
			{
			}
		};

	} //namespace curve

} //namespace dg
#endif //!CURVE_HB
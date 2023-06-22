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
			dg::math::Function X_;
			dg::math::Function Y_;
			dg::math::Function Z_;

		public:
			//double operator()(T var);
		};

	} //namespace curve

} //namespace dg
#endif //!CURVE_H

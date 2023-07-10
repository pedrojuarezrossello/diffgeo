#ifndef DG_CURVE_INTERFACE_H
#define DG_CURVE_INTERFACE_H
#include "../utils/real_function.h"
#include "../vector/vector.h"
#include <type_traits>

namespace dg {

	namespace curve {

		template <typename T>
		class RegularCurve;

		template<typename CurveImpl, typename T>
		class CurveInterface
		{
		protected:

			dg::math::Function<T,T> X_;
			dg::math::Function<T,T> Y_;
			dg::math::Function<T,T> Z_;

			template<int Order, typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			dg::vector::Vector<V> derivative_(V var) const
			{
				dg::vector::Vector<V> componentwiseDerivative(
					X_.derivative<Order>(var),
					Y_.derivative<Order>(var),
					Z_.derivative<Order>(var)
				);

				return componentwiseDerivative;
			}

			CurveInterface() {}

		public:

			//constructor

			CurveInterface(const dg::math::Function<T,T>& X, const dg::math::Function<T,T>& Y, const dg::math::Function<T,T>& Z)
				: X_(X), Y_(Y), Z_(Z) {}

			CurveInterface(const RegularCurve<T>& curve) : X_(std::move(curve.getX())), Y_(std::move(curve.getY())), Z_(std::move(curve.getZ())) {}
			
			//evaluation
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			dg::vector::Vector<V> operator()(V var)
			{
				dg::vector::Vector<V> evaluation(X_(var), Y_(var), Z_(var));
				return evaluation;
			}

			//derivative function
			template<int Order, typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			std::function<dg::vector::Vector<V>(V)> derivative() const
			{
				auto derivativeLambda = [this](V var) {
					return this->derivative_<Order>(var);
				};

				return derivativeLambda;
			}

			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			V curvature(V var) 
			{
				return static_cast<CurveImpl*>(this)->curvature_(var);
			}

			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			V torsion(V var)
			{
				return static_cast<CurveImpl*>(this)->torsion_(var);
			}

			template<typename V, 
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			V length(V a, V b)
			{
				return static_cast<CurveImpl*>(this)->length_(a,b);
			}

			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			dg::vector::Vector<V> unitTangent(V var)
			{
				return static_cast<CurveImpl*>(this)->unitTangent_(var);
			}

		};

	}
}
#endif //!DG_CURVE_INTERFACE_H
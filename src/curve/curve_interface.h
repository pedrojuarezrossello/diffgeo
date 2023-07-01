#ifndef CURVE_INTERFACE_H
#define CURVE_INTERFACE_H
#include <type_traits>
#include "../utils/real_function.h"
#include "../vector/vector.h"

namespace dg {

	namespace curve {

		template<typename CurveImpl, typename T>
		class CurveInterface
		{
		protected:

			dg::math::Function<T> X_;
			dg::math::Function<T> Y_;
			dg::math::Function<T> Z_;

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

			//unit tangent derivative function
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			auto unitTangentDerivative_(const dg::math::Function<T>& component, V var) const
			{
				auto normFunction = [this](auto x) {
					using std::sqrt;
					return sqrt(this->X_.derivative<1>(x) * this->X_.derivative<1>(x)
						+ this->Y_.derivative<1>(x) * this->Y_.derivative<1>(x)
						+ this->Z_.derivative<1>(x) * this->Z_.derivative<1>(x));
				};

				auto derivativeComponent = [&component](auto x) {
					return component.derivative<1>(x);
				};

				auto componentwiseUnitTangent = [&normFunction, &derivativeComponent](auto x)
				{
					return derivativeComponent(x) / normFunction(x);
				};

				auto const epsilonedValue = boost::math::differentiation::make_fvar<V, 1>(var);
				auto derivatives = componentwiseUnitTangent(epsilonedValue);
				return derivatives;

			}

			CurveInterface() {}

		public:

			//constructor

			CurveInterface(const dg::math::Function<T>& X, const dg::math::Function<T>& Y, const dg::math::Function<T>& Z)
				: X_(X), Y_(Y), Z_(Z) {}

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
			V unitTangent(V var)
			{
				return static_cast<CurveImpl*>(this)->unitTangent_(var);
			}

		};

	}
}
#endif //!CURVE_INTERFACE_H
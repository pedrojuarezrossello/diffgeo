#ifndef CURVE_H
#define CURVE_H
#include <type_traits>
#include "../utils/real_function.h"
#include "../vector/vector.h"
#include "../utils/operations.h"

namespace dg {

	namespace curve {

		template<typename T>
		class Curve
		{
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

		public:

			//constructor

			Curve() = default;

			Curve(const dg::math::Function<T>& X, const dg::math::Function<T>& Y, const dg::math::Function<T>& Z)
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

			//curvature (pointwise)
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			V curvature(V var) 
			{
				auto firstDerivative(derivative_<1>(var));
				auto crossWithSecondDerivative(dg::vector::cross_product(derivative_<2>(var), firstDerivative));
				auto firstDerivativeNorm{ firstDerivative.norm() };
				return crossWithSecondDerivative.norm() / (firstDerivativeNorm * firstDerivativeNorm * firstDerivativeNorm);

			}

			//torsion (pointwise)
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			V torsion(V var)
			{
				auto firstDerivative(derivative_<1>(var));
				auto secondDerivative(derivative_<2>(var));
				auto crossDerivatives(dg::vector::cross_product(secondDerivative,firstDerivative));
				auto crossNorm{ crossDerivatives.norm() };
				return (crossDerivatives * derivative_<3>(var)) / (crossNorm * crossNorm);
			}

		};

	} //namespace curve

} //namespace dg
#endif //!CURVE_HB
#ifndef CURVE_H
#define CURVE_H
#include <type_traits>
#include <boost/math/quadrature/gauss.hpp>
#include "../curve/curve_interface.h"
#include "../utils/real_function.h"
#include "../vector/vector.h"
#include "../utils/operations.h"
#include "../utils/maths.h"

namespace dg {

	namespace curve {

		template<typename T>
		class RegularCurve : public CurveInterface<RegularCurve<T>,T>
		{
		
		public:
			
			using CurveInterface<RegularCurve<T>, T>::CurveInterface;

			//curvature (pointwise)
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			V curvature_(V var) 
			{
				auto firstDerivative(this->derivative_<1>(var));
				auto crossWithSecondDerivative(dg::vector::cross_product(this->derivative_<2>(var), firstDerivative));
				auto firstDerivativeNorm{ firstDerivative.norm() };
				return crossWithSecondDerivative.norm() / (firstDerivativeNorm * firstDerivativeNorm * firstDerivativeNorm);

			}

			//torsion (pointwise)
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			V torsion_(V var)
			{
				auto firstDerivative(this->derivative_<1>(var));
				auto secondDerivative(this->derivative_<2>(var));
				auto crossDerivatives(dg::vector::cross_product(secondDerivative,firstDerivative));
				auto crossNorm{ crossDerivatives.norm() };
				return (crossDerivatives * this->derivative_<3>(var)) / (crossNorm * crossNorm);
			}

			//length
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			V length_(V a, V b)
			{
				auto normFunction = [this](auto x) {
					using std::sqrt;
					return sqrt(this->X_.derivative<1>(x) * this->X_.derivative<1>(x)
						+ this->Y_.derivative<1>(x) * this->Y_.derivative<1>(x)
						+ this->Z_.derivative<1>(x) * this->Z_.derivative<1>(x));
				};

				return boost::math::quadrature::gauss<V, 7>::integrate(normFunction, a, b);
			}

			//unit tangent vector
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			dg::vector::Vector<V> tangentVector_(V var)
			{
				dg::vector::Vector<V> firstDerivative(this->derivative_<1>(var));
				firstDerivative.normalise();
				return firstDerivative;
			}

			//reparametrise
			template<typename Functor>
			void reparametrise(Functor&& reparametrisation)
			{
				this->X_ = dg::math::compose<T>(this->X_, std::forward<Functor>(reparametrisation));
				this->Y_ = dg::math::compose<T>(this->Y_, std::forward<Functor>(reparametrisation));
				this->Z_ = dg::math::compose<T>(this->Z_, std::forward<Functor>(reparametrisation));
			}

			//total curvature
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			V totalCurvature(V startingPoint, V period) 
			{
				auto curvatureFunction = [startingPoint, period](auto x) {
					return curvature_(x);
				};

				return boost::math::quadrature::gauss<V, 7>::integrate(curvatureFunction, startingPoint, startingPoint+period);
			}

		};

	} //namespace curve

} //namespace dg
#endif //!CURVE_HB
#ifndef UNIT_CURVE_H
#define UNIT_CURVE_H
#include <type_traits>
#include "../curve/curve_interface.h"
#include "../curve/curve.h"
#include "../utils/real_function.h"
#include "../vector/vector.h"
#include "../utils/operations.h"
#include <iostream>

namespace dg {

	namespace curve {

		template<typename T>
		class UnitCurve : public CurveInterface<UnitCurve<T>,T>
		{

		public:

			using CurveInterface<UnitCurve<T>, T>::CurveInterface;

			//curvature (pointwise)
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			V curvature_(V var)
			{
				return std::sqrt(this->X_.derivative<2>(var) * this->X_.derivative<2>(var)
					+ this->Y_.derivative<2>(var) * this->Y_.derivative<2>(var)
					+ this->Z_.derivative<2>(var) * this->Z_.derivative<2>(var));
			}

			//torsion (pointwise)
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			V torsion_(V var)
			{
				auto firstDerivative(this->derivative_<1>(var));
				auto secondDerivative(this->derivative_<2>(var));
				auto crossDerivatives(dg::vector::cross_product(secondDerivative, firstDerivative));
				auto crossNorm{ crossDerivatives.norm() };
				return (crossDerivatives * this->derivative_<3>(var)) / (crossNorm * crossNorm);
			}

			//length
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			V length_(V a, V b) 
			{
				return b - a;
			
			}

			//unit tangent vector
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			dg::vector::Vector<V> tangentVector_(V var) 
			{
				dg::vector::Vector<V> firstDerivative(this->derivative_<1>(var));
				return firstDerivative;
			}

			//principal normal
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			dg::vector::Vector<V> principalNormal(V var) 
			{
				dg::vector::Vector<V> tangentDerivative(this->derivative_<2>(var));
				auto curv{ curvature_(var) };
				dg::vector::Vector<V> normal(curv * tangentDerivative);
				return normal;
			}

			
		};

	} //namespace curve

} //namespace dg


#endif // !UNIT_CURVE_H

/** @file unit_curve.h
 *	@ingroup Curves
 *  @brief Curve class for unit-speed curves.
 *
 *	It provides an implementation of the Curve interface
 *  for unit-speed curves.
 *
 *  @author Pedro Juarez Rossello (pedrojuarezrossello)
 *  @bug No known bugs.
 */

#ifndef UNIT_CURVE_H
#define UNIT_CURVE_H
#include <type_traits>
#include "../curve/curve_interface.h"
#include "../utils/real_function.h"
#include "../vector/vector.h"
#include "../vector/operations.h"
#include "../surface/plane.h"

namespace dg {

	namespace curve {

		template <typename T>
		class RegularCurve;

		/**
		* @headerfile unit_curve.h "src/curve/unit_curve.h"
		*
		* @tparam T The library has defined a type `dg::param<U>`
		* for U a flaoting point (underlying) type
		* which has to be passed as a template argument for the
		* component function of a curve to allow autodifferentiation.
		*/
		template<typename T>
		class UnitCurve : public CurveInterface<UnitCurve<T>,T>
		{
			//Declare interface friend class to hide implementation
			template<typename CurveImpl, typename U>
			friend class CurveInterface;

			//Curvature (pointwise) implementation
			template<typename V>
			V curvature_(V var)
			{
				return std::sqrt(this->X_.derivative<2>(var) * this->X_.derivative<2>(var)
					+ this->Y_.derivative<2>(var) * this->Y_.derivative<2>(var)
					+ this->Z_.derivative<2>(var) * this->Z_.derivative<2>(var));
			}

			//Torsion (pointwise) implementation
			template<typename V>
			V torsion_(V var)
			{
				auto firstDerivative(this->derivative_<1>(var));
				auto secondDerivative(this->derivative_<2>(var));
				auto crossDerivatives(dg::vector::cross_product(secondDerivative, firstDerivative));
				auto crossNorm{ crossDerivatives.norm() };
				return (crossDerivatives * this->derivative_<3>(var)) / (crossNorm * crossNorm);
			}

			//Length implementation
			template<typename V>
			V length_(V a, V b)
			{
				return b - a;

			}

			//Unit tangent vector (pointwise) implementation
			template<typename V>
			dg::vector::Vector<V> unitTangent_(V var)
			{
				dg::vector::Vector<V> firstDerivative(this->derivative_<1>(var));
				return firstDerivative;
			}

			//Principal normal util
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			dg::vector::Vector<V> principalNormal_(V var)
			{
				dg::vector::Vector<V> tangentDerivative(this->derivative_<2>(var));
				auto curv{ curvature_(var) };
				dg::vector::Vector<V> normal(curv * tangentDerivative);
				return normal;
			}

			//Binormal normal util
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			dg::vector::Vector<V> binormal_(V var)
			{
				dg::vector::Vector<V> binormalVec(dg::vector::cross_product(unitTangent_(var), principalNormal(var)));
				return binormalVec;
			}

			//Osculating plane util
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			dg::surf::Plane<V, V> osculatingPlane_(V var)
			{
				dg::surf::Plane<V, V> osculating(unitTangent_(var), principalNormal(var), this->operator()(var));
				return osculating;
			}

			//Rectifying plane util
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			dg::surf::Plane<V, V> rectifyingPlane_(V var)
			{
				dg::surf::Plane<V, V> rectifying(binormal(var), unitTangent_(var), this->operator()(var));
				return rectifying;
			}

			//Normal plane util
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			dg::surf::Plane<V, V> normalPlane_(V var)
			{
				dg::surf::Plane<V, V> normal(principalNormal(var), binormal(var), this->operator()(var));
				return normal;
			}

		public:

			using CurveInterface<UnitCurve<T>, T>::CurveInterface;

			//Constructor

			UnitCurve(const RegularCurve<T>& curve) : CurveInterface<UnitCurve<T>, T>(curve) {}

			/**
			* @brief Principal normal at parameter value `var`.
			*
			* @tparam V Floating point type of `var`.
			* @param var Floating point value at which we compute the principal normal vector.
			*
			* @return Vector of type V.
			*/
			template<typename V>
			dg::vector::Vector<V> principalNormal(V var) 
			{
				return principalNormal_(var);
			}

			/**
			* @brief Binormal at parameter value `var`.
			*
			* @tparam V Floating point type of `var`.
			* @param var Floating point value at which we compute the binormal vector.
			*
			* @return Vector of type V.
			*/
			template<typename V>
			dg::vector::Vector<V> binormal(V var)
			{
				return binormal_(var);
			}

			/**
			* @brief Osculating plane at parameter value `var`.
			*
			* @details It's spanned by the unit tangent and the principal normal.
			* 
			* @tparam V Floating point type of `var`.
			* @param var Floating point value at which we compute the osculating plane.
			*
			* @return Plane object.
			*/
			template<typename V>
			dg::surf::Plane<V, V> osculatingPlane(V var)
			{
				return osculatingPlane_(var);
			}
			
			/**
			* @brief Rectifying plane at parameter value `var`.
			*
			* @details It's spanned by the binormal and the unit tangent.
			* 
			* @tparam V Floating point type of `var`.
			* @param var Floating point value at which we compute the rectifying plane.
			*
			* @return Plane object.
			*/
			template<typename V>
			dg::surf::Plane<V, V> rectifyingPlane(V var)
			{
				return rectifyingPlane_(var);
			}

			/**
			* @brief Normal plane at parameter value `var`.
			*
			* @details It's spanned by the binormal and the principal normal.
			* 
			* @tparam V Floating point type of `var`.
			* @param var Floating point value at which we compute the normal plane.
			*
			* @return Plane object.
			*/
			template<typename V>
			dg::surf::Plane<V, V> normalPlane(V var)
			{
				return normalPlane_(var);
			}
			
		};

	} //namespace curve

} //namespace dg


#endif // !UNIT_CURVE_H

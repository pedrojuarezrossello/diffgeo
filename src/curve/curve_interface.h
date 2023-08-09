 /** @file curve_interface.h
 *	@ingroup Curves
 *  @brief Interface class for curves.
 *
 *	It provides the basic interface for the curve
 *	classes RegularCurve and UnitCurve. It uses the
 *	CRTP to achieve a sort of reversed (static) polymorphism,
 *	so the functions that one must call from client code are
 *	specified here. 
 * 
 *  @author Pedro Juarez Rossello (pedrojuarezrossello)
 *  @bug No known bugs.
 */

#ifndef DG_CURVE_INTERFACE_H
#define DG_CURVE_INTERFACE_H
#include "../utils/real_function.h"
#include "../vector/vector.h"
#include <type_traits>
#include <iostream>

namespace dg {

	namespace curve {

		template <typename T>
		class RegularCurve;

		/**
		* @headerfile curve_interface.h "src/curve/curve_interface.h"
		* 
		*/
		template<typename CurveImpl, typename T>
		class CurveInterface
		{
		
			//Evaluation util method
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			dg::vector::Vector<V> evaluate_(V var)
			{
				dg::vector::Vector<V> evaluation(X_(var), Y_(var), Z_(var));
				return evaluation;
			}

			//Curvature util method
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			V curvatureUtil_(V var)
			{
				return static_cast<CurveImpl*>(this)->curvature_(var);
			}

			//Torsion util method
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			V torsionUtil_(V var)
			{
				return static_cast<CurveImpl*>(this)->torsion_(var);
			}

			//Length util method
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			V lengthUtil_(V a, V b)
			{
				return static_cast<CurveImpl*>(this)->length_(a, b);
			}

			//Unit tangent util method
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			dg::vector::Vector<V> unitTangentUtil_(V var)
			{
				return static_cast<CurveImpl*>(this)->unitTangent_(var);
			}

		protected:

			/**
			* @brief X component function.
			* 
			* @tparam T The library has defined a type `dg::param<U>` 
			* for U a flaoting point (underlying) type
			* which has to be passed as a template argument for the 
			* component function of a curve to allow autodifferentiation.
			*/
			dg::math::Function<T,T> X_;

			/**
			* @brief X component function.
			*
			* @tparam T The library has defined a type `dg::param<U>`
			* for U a flaoting point (underlying) type
			* which has to be passed as a template argument for the
			* component function of a curve to allow autodifferentiation.
			*/
			dg::math::Function<T,T> Y_;

			/**
			* @brief X component function.
			*
			* @tparam T The library has defined a type `dg::param<U>`
			* for U a flaoting point (underlying) type
			* which has to be passed as a template argument for the
			* component function of a curve to allow autodifferentiation.
			*/
			dg::math::Function<T,T> Z_;
			
			//Derivative util method 
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

			//We hide the constructor 
			CurveInterface() {}

		public:

			//Constructors

			CurveInterface(const dg::math::Function<T,T>& X, const dg::math::Function<T,T>& Y, const dg::math::Function<T,T>& Z)
				: X_(X), Y_(Y), Z_(Z) {}

			CurveInterface(const RegularCurve<T>& curve) : X_(std::move(curve.getX())), Y_(std::move(curve.getY())), Z_(std::move(curve.getZ())) {}
			
			/**
			* @brief Evalutation at parameter value `var`.
			* 
			* @details It returns the result as a Vector of type V.
			* 
			* @tparam V Floating point type of the return Vector.
			* @param var Floating point value at which we evaluate the curve.
			* 
			* @return Vector of type V.
			*/
			template<typename V>
			dg::vector::Vector<V> operator()(V var)
			{
				return evaluate_<V>(var);
			}

			/**
			* @brief Derivative function of order Order.
			*
			* @details It returns an `std::function` object
			* that outputs the derivative as a Vector of type V 
			* for a given value of type V.
			* Note Order is a template parameter.
			*
			* @tparam Order Order (int) of the derivative.
			* @tparam V Floating point type of the return Vector.
			* 
			* @return `std::function` containing the derivative function
			*/
			template<int Order, typename V>
			std::function<dg::vector::Vector<V>(V)> derivative() const
			{
				auto derivativeLambda = [this](V var) {
					return this->derivative_<Order>(var);
				};

				return derivativeLambda;
			}

			/**
			* @brief Curvature at parameter value `var`.
			*
			* @tparam V Floating point type of `var`.
			* @param var Floating point value at which we compute the curvature.
			* 
			* @return Curvature as a floating point type V.
			*/
			template<typename V>
			V curvature(V var) 
			{
				return curvatureUtil_<V>(var);
			}

			/**
			* @brief Torsion at parameter value `var`.
			*
			* @tparam V Floating point type of `var`.
			* @param var Floating point value at which we compute the torsion.
			*
			* @return Torsion as a floating point type V.
			*/
			template<typename V>
			V torsion(V var)
			{
				return torsionUtil_<V>(var);
			}

			/**
			* @brief Length between `a` and `b`.
			*
			* @details It uses `boost::math::quadrature::gauss` 
			* to integrate the length function.
			*
			* @tparam V Floating point type of `a` and `b`.
			* @param a Start of range of integration.
			* @param b End of range of integration.
			*
			* @return Torsion as a floating point type V.
			*/
			template<typename V>
			V length(V a, V b)
			{
				return lengthUtil_<V>(a,b);
			}

			/**
			* @brief Unit tangent vector at parameter value `var`.
			*
			* @tparam V Floating point type of `var`.
			* @param var Floating point value at which we compute the unit tangent vector.
			*
			* @return Vector of type V.
			*/
			template<typename V>
			dg::vector::Vector<V> unitTangent(V var)
			{
				return unitTangentUtil_<V>(var);
			}

		};

	}
}
#endif //!DG_CURVE_INTERFACE_H
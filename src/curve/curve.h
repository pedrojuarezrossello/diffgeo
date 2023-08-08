/** @file curve.h
 *	@ingroup Curves
 *  @brief Curve class for non unit-length curves.
 *
 *	It provides an implementation of the Curve interface
 *  for curves that aren't necessarily unit-length.
 *
 *  @author Pedro Juarez Rossello (pedrojuarezrossello)
 *  @bug No known bugs.
 */

#ifndef DG_CURVE_H
#define DG_CURVE_H
#include "../curve/curve_interface.h"
#include "../utils/real_function.h"
#include "../vector/vector.h"
#include "../vector/operations.h"
#include "../utils/maths.h"
#include "../curve/unit_curve.h"
#include <boost/math/quadrature/gauss.hpp>
#include <type_traits>
#include <memory>
#include <iostream>

namespace dg {

	namespace curve {

		/**
		* @headerfile curve_interface.h "src/curve/curve_interface.h"
		* 
		* @tparam T The library has defined a type `dg::param<U>` 
		* for U a flaoting point (underlying) type
		* which has to be passed as a template argument for the 
		* component function of a curve to allow autodifferentiation.
		*/
		template<typename T>
		class RegularCurve : public CurveInterface<RegularCurve<T>,T>
		{
			//Declare interface friend class to hide implementation
			template<typename CurveImpl,typename U>
			friend class CurveInterface;

			//Curvature (pointwise) implementation
			template<typename V>
			V curvature_(V var)
			{
				auto firstDerivative(this->derivative_<1>(var));
				auto crossWithSecondDerivative(dg::vector::cross_product(this->derivative_<2>(var), firstDerivative));
				auto firstDerivativeNorm{ firstDerivative.norm() };
				return crossWithSecondDerivative.norm() / (firstDerivativeNorm * firstDerivativeNorm * firstDerivativeNorm);
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
				auto normFunction = [this](auto x) {
					using std::sqrt;
					return sqrt(this->X_.derivative<1>(x) * this->X_.derivative<1>(x)
						+ this->Y_.derivative<1>(x) * this->Y_.derivative<1>(x)
						+ this->Z_.derivative<1>(x) * this->Z_.derivative<1>(x));
				};

				return boost::math::quadrature::gauss<V, 7>::integrate(normFunction, a, b);
			}

			//Unit tangent vector (pointwise) implementation
			template<typename V>
			dg::vector::Vector<V> unitTangent_(V var)
			{
				dg::vector::Vector<V> firstDerivative(this->derivative_<1>(var));
				firstDerivative.normalise();
				return firstDerivative;
			}

			//Total curvature util
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			V totalCurvature_(V startingPoint, V period)
			{
				auto curvatureFunction = [startingPoint, period, this](auto x) {
					return this->curvature_(x);
				};

				return boost::math::quadrature::gauss<V, 7>::integrate(curvatureFunction, startingPoint, startingPoint + period);
			}

		public:
			
			using CurveInterface<RegularCurve<T>, T>::CurveInterface;

			/**
			* @brief Getter function for the first component.
			* 
			* @return `const`-reference to the function.
			*/
			const auto& getX() const { return this->X_; }

			/**
			* @brief Getter function for the second component.
			* 
			* @return `const`-reference to the function.
			*/
			const auto& getY() const { return this->Y_; }

			/**
			* @brief Getter function for the third component.
			* 
			* @return `const`-reference to the function.
			*/
			const auto& getZ() const { return this->Z_; }
			
			/**
			* @brief Reparametrise the curve.
			* 
			* @details The function argument has to be a template function
			* that's passed a template parameter `dg::param<V>`
			* for the desired underlying floating point type.
			* For instance, if we wanted to use the reparametrisation x/2
			* then we'd have to define ``template<typename T> auto repam(const T& x) {return 0.5*x;}``,
			* and then pass repam<dg::param_d<>> to the method.
			* 
			* @tparam Functor Function type of the method parameter.
			* @param reparametrisation Function representing the desired reparametrisation.
			* It needs to follow the rules above.
			* 
			* @return void.
			*/
			template<typename Functor>
			void reparametrise(Functor&& reparametrisation)
			{
				this->X_ = dg::math::compose<Functor,T,T>(this->X_, std::forward<Functor>(reparametrisation));
				this->Y_ = dg::math::compose<Functor,T,T>(this->Y_, std::forward<Functor>(reparametrisation));
				this->Z_ = dg::math::compose<Functor,T,T>(this->Z_, std::forward<Functor>(reparametrisation));
			}


			/**
			* @brief Reparametrise the curve to unit speed.
			*
			* @details Note that the unit speed parametrisation has to be provided 
			* as an argument. The function argument has to be a template function
			* that's passed a template parameter `dg::param<V>`
			* for the desired underlying floating point type.
			* For instance, if we wanted to use the reparametrisation x/2
			* then we'd have to define ``template<typename T> auto repam(const T& x) {return 0.5*x;}``,
			* and then pass repam<dg::param_d<>> to the method.
			*
			* @tparam Functor Function type of the method parameter.
			* @param reparametrisation Function representing the desired reparametrisation.
			* It needs to follow the rules above.
			*
			* @return An `std::unique_ptr` to a new UnitCurve object.
			*/
			template<typename Functor>
			std::unique_ptr<UnitCurve<T>> unitSpeedParametrisation(Functor&& param)
			{
				reparametrise(std::forward<Functor>(param));
				return std::make_unique<UnitCurve<T>>(*this);
			}

			/**
			* @brief Total curvature of a closed curve
			* 
			* @details Note that this value is only meaningful for closed curves.
			* 
			* @tparam V Floating point type we'd like the result to be.
			* @param startingPoint Initial value (can be anything so long as 
			* the curve is well-defined at startingPoint + period).
			* @param period Period of the closed curve.
			* 
			* @return V Floating point type.
			*/
			template<typename V>
			V totalCurvature(V startingPoint, V period) 
			{
				return totalCurvature_(startingPoint, period);
			}
		};

	} //namespace curve

} //namespace dg
#endif //!DG_CURVE_H
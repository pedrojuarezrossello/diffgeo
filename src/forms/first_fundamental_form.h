/** @file first_fundamental_form.h
 *	@ingroup Surfaces
 *  @brief First fundamental form class.
 *
 *  @author Pedro Juarez Rossello (pedrojuarezrossello)
 *  @bug No known bugs.
 */

#ifndef DG_FIRST_FUNDAMENTAL_FORM_H
#define DG_FIRST_FUNDAMENTAL_FORM_H
#include "../forms/two_form_interface.h"
#include "../surface/surface.h"
#include "../utils/real_function.h"

namespace dg {

	namespace form {

		/**
		* @headerfile first_fundamental_form.h "src/forms/first_fundamental_form.h" 
		* 
		* @tparam T Floating point type of component functions.
		*/
		template<typename T>
		class FirstFundamentalForm : public TwoFormInterface<T>
		{
		
			//Line element util method
			auto computeLineElement_(T X_comp, T Y_comp, T X_der, T Y_der) const
			{
				return this->A_(X_comp, Y_comp) * X_der * X_der
					+ 2.0 * this->B_(X_comp, Y_comp) * X_der * Y_der
					+ this->C_(X_comp, Y_comp) * Y_der * Y_der;
			}

		public:

			//Constructors
			FirstFundamentalForm() = default;

			FirstFundamentalForm( std::function<T(T, T)> A,  std::function<T(T, T)> B,  std::function<T(T, T)> C)
				: TwoFormInterface<T>(A,B,C) {}

			/**
			* @brief Line element function
			* 
			* @details It returns the coefficient f(t) in the line element ds = f(t) dt
			* for a given curve on the surface we're dealing with.
			* 
			* @tparam CurveParam Library-specific type `dg::param<U>` 
			* for U a floating point (underlying) type
			* which has to be passed as a template argument for the 
			* component function of a curve to allow autodifferentiation.
			* See example paraboloid_curve_length.cpp.
			* Also note that the curve only has two components as it lives 
			* on the surface and not on the ambient 3D space.
			* 
			* @param X_component First component of the curve.
			* @param Y_component Second component of the curve.
			* 
			* @return A lambda returning the coefficient of the the line element at each value.
			*/
			template<typename CurveParam>
			auto lineElement(const dg::math::Function<CurveParam, CurveParam>& X_component,
				const dg::math::Function<CurveParam, CurveParam>& Y_component) const
			{
				auto lineElem = [this, X_component, Y_component](auto x) {
					using std::sqrt;
					return sqrt(this->computeLineElement_(X_component.derivative<0>(x), Y_component.derivative<0>(x), X_component.derivative<1>(x), Y_component.derivative<1>(x)));
				};

				return lineElem;
			}

			/**
			* @brief Area element function
			*
			* @details It returns the coefficient f(u,v) 
			* in the area element dA = f(u,v) dudv.
			*
			* @return A lambda returning the coefficient 
			* of the area element at each pair of values u and v (of type T,
			* same as the form.
			*/
			std::function<T(T,T)> areaElement() const
			{
				auto areaElement = [this](T var1, T var2) -> T {
					using std::sqrt;
					return sqrt(
						this->A_(var1, var2) * this->C_(var1, var2) - this->B_(var1, var2) * this->B_(var1, var2)
					);
				};

				return areaElement;
			}

		};

	} //namespace form

} //namespace dg

#endif // !FIRST_FUNDAMENTAL_FORM_H

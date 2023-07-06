#ifndef FIRST_FUNDAMENTAL_FORM_H
#define FIRST_FUNDAMENTAL_FORM_H
#include <iostream>
#include "../forms/two_form_interface.h"
#include "../surface/surface.h"
#include "../utils/real_function.h"

namespace dg {

	namespace form {

		template<typename T>
		class FirstFundamentalForm : public TwoFormInterface<T>
		{
		
			auto computeLineElement_(T X_comp, T Y_comp, T X_der, T Y_der) const
			{
				return this->A_(X_comp, Y_comp) * X_der * X_der
					+ 2.0 * this->B_(X_comp, Y_comp) * X_der * Y_der
					+ this->C_(X_comp, Y_comp) * Y_der * Y_der;
			}


		public:

			FirstFundamentalForm( std::function<T(T, T)> A,  std::function<T(T, T)> B,  std::function<T(T, T)> C)
				: TwoFormInterface<T>(A,B,C) {}

			//line element
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

			//area element
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

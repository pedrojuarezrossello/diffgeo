#ifndef FIRST_FUNDAMENTAL_FORM_H
#define FIRST_FUNDAMENTAL_FORM_H
#include <iostream>
#include "../forms/two_form_interface.h"
#include "../surface/surface.h"
#include "../utils/real_function.h"

namespace dg {

	namespace form {

		template<typename T>
		class FirstFundamentalForm 
		{
			std::function<T(T, T)> A_;
			std::function<T(T, T)> B_;
			std::function<T(T, T)> C_;
		
			auto computeLineElement_(T X_comp, T Y_comp, T X_der, T Y_der) const
			{
				return A_(X_comp, Y_comp) * X_der * X_der
					+ 2.0 * B_(X_comp, Y_comp) * X_der * Y_der
					+ C_(X_comp, Y_comp) * Y_der * Y_der;
			}


		public:

			FirstFundamentalForm( std::function<T(T, T)> A,  std::function<T(T, T)> B,  std::function<T(T, T)> C)
				: A_(A), B_(B), C_(C) {}

			/*template<typename _Ret, typename... Args>
			FirstFundamentalForm(const dg::surf::Surface<_Ret, Args...> surface) : A_(nullptr), 
				B_(nullptr),
				C_(nullptr) {
				A_ = surface.getE_<T>();
				B_ = surface.getF_<T>();
				C_ = surface.getG_<T>();
				std::cout << "A is " << A_(1.9, 0.3) << std::endl;
			}*/

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
						A_(var1, var2) * C_(var1, var2) - B_(var1, var2) * B_(var1, var2)
					);
				};

				return areaElement;
			}

		};

	} //namespace form

} //namespace dg

#endif // !FIRST_FUNDAMENTAL_FORM_H

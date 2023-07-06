#ifndef SECOND_FUNDAMENTAL_FORM_H
#define SECOND_FUNDAMENTAL_FORM_H
#include "../forms/two_form_interface.h"

namespace dg {

	namespace form {

		template<typename T>
		class SecondFundamentalForm : public TwoFormInterface<T>
		{

		public:

			template<typename _Ret, typename... Args>
			SecondFundamentalForm(const dg::surf::Surface<_Ret, Args...> surface) : TwoFormInterface(std::move(surface.getL_<T, T>()),
				std::move(surface.getM_<T, T>()),
				std::move(surface.getN_<T, T>())) {}

		};

	} //namespace form

} //namespace dg

#endif // !SECOND_FUNDAMENTAL_FORM_H

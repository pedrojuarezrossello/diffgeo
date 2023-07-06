#ifndef SECOND_FUNDAMENTAL_FORM_H
#define SECOND_FUNDAMENTAL_FORM_H
#include "../forms/two_form_interface.h"

namespace dg {

	namespace form {

		template<typename T>
		class SecondFundamentalForm : public TwoFormInterface<T>
		{

		public:

			SecondFundamentalForm(std::function<T(T, T)> A, std::function<T(T, T)> B, std::function<T(T, T)> C)
				: TwoFormInterface<T>(A, B, C) {}

		};

	} //namespace form

} //namespace dg

#endif // !SECOND_FUNDAMENTAL_FORM_H

/** @file second_fundamental_form.h
 *	@ingroup Surfaces
 *  @brief Second fundamental form class.
 *
 *  @author Pedro Juarez Rossello (pedrojuarezrossello)
 *  @bug No known bugs.
 */

#ifndef DG_SECOND_FUNDAMENTAL_FORM_H
#define DG_SECOND_FUNDAMENTAL_FORM_H
#include "../forms/two_form_interface.h"

namespace dg {

	namespace form {

		/**
		* @headerfile second_fundamental_form.h "src/forms/second_fundamental_form.h"
		*
		* @tparam T Floating point type of component functions.
		*/
		template<typename T>
		class SecondFundamentalForm : public TwoFormInterface<T>
		{

		public:

			//Constructor 

			SecondFundamentalForm(std::function<T(T, T)> A, std::function<T(T, T)> B, std::function<T(T, T)> C)
				: TwoFormInterface<T>(A, B, C) {}

		};

	} //namespace form

} //namespace dg

#endif // !DG_SECOND_FUNDAMENTAL_FORM_H

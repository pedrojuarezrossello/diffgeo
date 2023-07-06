#ifndef TWO_FORM_INTERFACE_H
#define TWO_FORM_INTERFACE_H
#include "../utils/real_function.h"

namespace dg {

	namespace form {

		template<typename T>
		class TwoFormInterface 
		{

		protected:

			std::function<T(T, T)> A_;
			std::function<T( T, T)> B_;
			std::function<T( T, T)> C_;
			int num;
			
			template<typename T>
			auto computeLineElement_(T X_comp, T Y_comp, T X_der, T Y_der) const
			{
				return A_(X_comp, Y_comp) * X_der * X_der
					+ 2.0 * B_(X_comp, Y_comp) * X_der * Y_der
					+ C_(X_comp, Y_comp) * Y_der * Y_der;
			}

			TwoFormInterface() {} //hide interface constructor

		public:

			TwoFormInterface(const std::function<T(T, T)>& A, const std::function<T(T, T)>& B, const std::function<T(T, T)>& C, int a)
				: A_(A), B_(B), C_(C), num(a) {}

			
			
		};

	} //namespace form

} //namespace dg

#endif // !TWO_FORM_INTERFACE_H


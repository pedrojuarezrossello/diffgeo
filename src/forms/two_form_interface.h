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
		

			TwoFormInterface() {} //hide interface constructor

		public:

			TwoFormInterface(std::function<T(T, T)> A, std::function<T(T, T)> B, std::function<T(T, T)> C)
				: A_(A), B_(B), C_(C) {}
			
		};

	} //namespace form

} //namespace dg

#endif // !TWO_FORM_INTERFACE_H


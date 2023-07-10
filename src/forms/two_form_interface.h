#ifndef DG_TWO_FORM_INTERFACE_H
#define DG_TWO_FORM_INTERFACE_H
#include "../utils/real_function.h"
#include "../forms/vector_field.h"

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

			template<typename U>
			auto innerProduct(const VectorField<U>& vector1, const VectorField<U>& vector2, U var1, U var2)
			{
				return A_(var1, var2) * vector1.compute_U(var1, var2) * vector2.compute_U(var1, var2)
					+ B_(var1, var2) * (vector1.compute_U(var1, var2) * vector2.compute_V(var1, var2) + vector1.compute_V(var1, var2) * vector2.compute_U(var1, var2))
					+ C_(var1, var2) * vector1.compute_V(var1, var2) * vector2.compute_V(var1, var2);
			}
			
		};

	} //namespace form

} //namespace dg

#endif // !DG_TWO_FORM_INTERFACE_H


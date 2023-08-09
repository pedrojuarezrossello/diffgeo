/** @file two_form_interface.h
 *	@ingroup Surfaces
 *  @brief Interface class for 2-forms.
 *
 *	It provides the basic interface for the symmetric 2-form
 *	classes FirstFundamentalForm and SecondFundamentalForm. 
 *	
 *  @author Pedro Juarez Rossello (pedrojuarezrossello)
 *  @bug No known bugs.
 */

#ifndef DG_TWO_FORM_INTERFACE_H
#define DG_TWO_FORM_INTERFACE_H
#include "../utils/real_function.h"
#include "../forms/vector_field.h"

namespace dg {

	namespace form {

		/**
		* @headerfile two_form_interface.h "src/forms/two_form_interface.h"
		*/
		template<typename T>
		class TwoFormInterface 
		{

		protected:

			/**
			* @brief First component of the two-form
			* 
			* @details It's implemented as an `std::function` type,
			* so it can accept any callable object.
			*/
			std::function<T(T, T)> A_;

			/**
			* @brief Second component of the two-form
			*
			* @details It's implemented as an `std::function` type,
			* so it can accept any callable object.
			*/
			std::function<T(T, T)> B_;

			/**
			* @brief Third component of the two-form
			*
			* @details It's implemented as an `std::function` type,
			* so it can accept any callable object.
			*/
			std::function<T(T, T)> C_;
		
			//Hide constructors

			TwoFormInterface() {}

			TwoFormInterface(std::function<T(T, T)> A, std::function<T(T, T)> B, std::function<T(T, T)> C)
				: A_(A), B_(B), C_(C) {}

		public:

			/**
			* @brief Compute the inner product of two VectorField objects
			* 
			* @details It works out the inner product induced by the 2-form
			* of two vector fields at given parameter values `var1` and `var2`.
			* 
			* @tparam U Floating point underlying type.
			* @param vector1 First VectorField<U> object.
			* @param vector2 Second VectorField<U> object.
			* @param var1 First parameter value.
			* @param var2 Second parameter value.
			* 
			* @return A scalar of type U.
			*/
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


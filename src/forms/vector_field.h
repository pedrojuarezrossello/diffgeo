/** @file vector_field.h
 *	@ingroup Surfaces
 *  @brief Vector field class.
 * 
 *	@details It implements a vector field of the form 
 *	f(x,y) D/DX + g(x,y) D/DY, where f,g are smooth functions.
 *
 *  @author Pedro Juarez Rossello (pedrojuarezrossello)
 *  @bug No known bugs.
 */

#ifndef DG_VECTOR_FIELD_H
#define DG_VECTOR_FIELD_H
#include <functional>

namespace dg {

	namespace form {

		/**
		* @headerfile vector_field.h "src/forms/vector_field.h"
		* 
		* @tparam T Floating point type of component functions.
		*/
		template<typename T>
		class VectorField
		{
			//f(x,y) as above
			std::function<T(T, T)> U_component;

			//g(x,y) as above
			std::function<T(T, T)> V_component;

		public:

			/**
			* @brief Default constructor
			*/
			VectorField() = default;

			/**
			* @brief Constructor
			* 
			* @details It uses two function arguments to construct the data members.
			* 
			* @tparam Functor1 Callable object type 
			* @tparam Functor2 Callable object type 
			* @param U Function object of D/DX component
			* @param V Function object of D/DY component
			*/
			template<typename Functor1, typename Functor2>
			VectorField(Functor1&& U, Functor2&& V) : U_component(U), V_component(V) {}

			/**
			* @brief Evaluate first component function.
			* 
			* @return Scalar of type T.
			*/
			auto compute_U(T var1, T var2) const { return U_component(var1,var2); }

			/**
			* @brief Evaluae second component function.
			*
			* @return Scalar of type T.
			*/
			auto compute_V(T var1, T var2) const { return V_component(var1,var2); }

		};

	} //namespace form

} //namespace dg

#endif // !DG_VECTOR_FIELD_H

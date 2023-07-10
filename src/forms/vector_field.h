#ifndef VECTOR_FIELD_H
#define VECTOR_FIELD_H
#include <functional>

namespace dg {

	namespace form {

		template<typename T>
		class VectorField
		{
			std::function<T(T, T)> U_component;
			std::function<T(T, T)> V_component;

		public:

			VectorField() {}

			template<typename Functor1, typename Functor2>
			VectorField(Functor1&& U, Functor2&& V) : U_component(U), V_component(V) {}

			T compute_U(T var1, T var2) const { return U_component(var1, var2); }

			T compute_V(T var1, T var2) const { return V_component(var1, var2); }

		};

	} //namespace form

} //namespace dg

#endif // !VECTOR_FIELD_H

#ifndef TYPE_TRAITS_H
#define TYPE_TRAITS_H
#include <type_traits>
#include <boost/math/differentiation/autodiff.hpp>

namespace dg {

	namespace vector {

		template <typename CallableObject, typename... Args>
		class VectorExpression;

		template <typename T>
		class Vector;

		template <typename T>
		struct is_vector_or_expression
		{
			static constexpr bool value = false;
		};

		template<typename T>
		struct is_vector_or_expression<Vector<T>>
		{
			static constexpr bool value = true;
		};

		template<typename CallableObject, typename ...Args>
		struct is_vector_or_expression <VectorExpression<CallableObject, Args...>>
		{
			static constexpr bool value = true;
		};

		template<typename T>
		constexpr bool is_vector_or_expression_t = is_vector_or_expression < std::remove_cv_t < std::remove_reference_t<T>>>::value;

		template<typename T>
		struct is_vector
		{
			static constexpr bool value = false;
		};

		template<typename T>
		struct is_vector<Vector<T>>
		{
			static constexpr bool value = true;
		};

		template<typename T>
		constexpr bool is_vector_v = is_vector<std::remove_cv_t < std::remove_reference_t<T>>>::value;


	} //namespace vector

	namespace math {

		template<typename T,
			std::enable_if_t<std::is_floating_point_v<T>, std::nullptr_t> = nullptr>
		struct number_helper {
			using type = boost::math::differentiation::autodiff_v1::autodiff_fvar<T, 3>;
		};

		template< typename T>
		using num = typename number_helper<T>::type;

		template <typename T>
		struct is_autodiff_compatible
		{
			static constexpr bool value = false;
		};

		/*template<typename T>
		struct is_autodiff_compatible<typename num<T>>
		{
			static constexpr bool value = true;
		};

		template<typename U>
		constexpr bool is_autodiff_compatible_v = is_autodiff_compatible<std::remove_cv_t < std::remove_reference_t<U>>>::value;*/
	}

} //namespace dg


#endif //!TYPE_TRAITS_H
#ifndef TYPE_TRAITS_H
#define TYPE_TRAITS_H
#include <type_traits>


template <typename CallableObject, typename... Args>
class VectorExpression;

template <typename T>
class Vector;

template<class T>
struct remove_cvref
{
	typedef std::remove_cv_t<std::remove_reference_t<T>> type;
};

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

template<typename T, typename U,
	std::enable_if_t<std::is_floating_point_v<T>, std::nullptr_t> = nullptr>
struct reduce_to_float
{ ///????
};

#endif //!TYPE_TRAITS_H
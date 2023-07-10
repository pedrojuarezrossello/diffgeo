#ifndef TYPE_TRAITS_H
#define TYPE_TRAITS_H
#include <boost/math/differentiation/autodiff.hpp>
#include <type_traits>


namespace dg {

	namespace surf {

		template<typename _Ret, typename... Args> //forward declaration
		class Surface;

	} //namespace surf
	 
	namespace curve {

		template<typename T>
		class RegularCurve; //forward declaration

		template<typename T>
		class UnitCurve; //forward declaration

	} //namespace curve

	//Type traits compatible with Autodifferentiation BOOST LIB

	/*
	Parameter to pass in template instantiation 
	of Curve component functions, so if for example template<typename T> auto function(T x) {}
	is the X component of a curve, we'd pass into the Function constructor as function<param_d>
	*/
	template<typename T>
	struct curv_param_ {
		using type = boost::math::differentiation::autodiff_v1::autodiff_fvar<T, 3>;
	};

	template<typename T>
	using param = typename curv_param_<T>::type;

	template<typename T=void>
	using param_d = param<double>;

	template<typename T=void>
	using param_f = param<float>;

	template<typename T=void>
	using param_ld = param<long double>;

	/*
	Parameter to pass in template instantiation
	of Surface component functions, so if for example template<typename T, typename S> auto function(T x, S y) {}
	is the X component of a curve, we'd pass into the Function constructor as function<param1_d, param2_d>
	*/

	template<typename T>
	struct surf_1_param_ {
		using type = boost::math::differentiation::autodiff_v1::detail::fvar<T, 2>;
	};

	template<typename T>
	using param1 = typename surf_1_param_<T>::type;

	template<typename T=void>
	using param1_d = param1<double>;

	template<typename T=void>
	using param1_f = param1<float>;

	template<typename T=void>
	using param1_ld = param1<long double>;

	template<typename T>
	struct surf_2_param_ {
		using type = boost::math::differentiation::autodiff_v1::detail::fvar<boost::math::differentiation::autodiff_v1::detail::fvar<T, 2>, 0>;
	};

	template<typename T>
	using param2 = typename surf_2_param_<T>::type;

	template<typename T=void>
	using param2_d = param2<double>;

	template<typename T=void>
	using param2_f = param2<float>;

	template<typename T=void>
	using param2_ld = param2<long double>;

	/*
	* Component function for curves
	*/

	template<typename T>
	struct Component_function_curve_ {
		using type = dg::math::Function< boost::math::differentiation::autodiff_v1::detail::fvar<T, 3>, boost::math::differentiation::autodiff_v1::detail::fvar<T, 3>>;
	};

	template<typename T>
	using CurveComponent = typename Component_function_curve_<T>::type;

	template<typename T=void>
	using CurveComponent_d = CurveComponent<double>;

	template<typename T=void>
	using CurveComponent_f = CurveComponent<float>;

	template<typename T=void>
	using CurveComponent_ld = CurveComponent<long double > ;

	/*
	* Component function for surfaces
	*/

	template<typename T>
	struct Surf_param_ {
		using type = boost::math::differentiation::autodiff_v1::autodiff_fvar<T, 2>;
	};

	template< typename T>
	using surf_param__ = typename Surf_param_<T>::type;

	template<typename T>
	using surf_param_return_ = surf_param__<surf_param__<T>>;

	template<typename T>
	struct Component_function_surface_ {
		using type = dg::math::Function<surf_param_return_<T>, param1<T>, param2<T>>;
	};

	template<typename T>
	using SurfaceComponent = typename Component_function_surface_<T>::type;

	template<typename T=void>
	using SurfaceComponent_d = SurfaceComponent<double>;

	template<typename T=void>
	using SurfaceComponent_f = SurfaceComponent<float>;

	template<typename T=void>
	using SurfaceComponent_ld = SurfaceComponent<long double>;

	/*
	* Surfaces
	*/

	template<typename T>
	struct Surface_ {
		using type = dg::surf::Surface<surf_param_return_<T>, param1<T>, param2<T>>;
	};

	template<typename T>
	using Surface = typename Surface_<T>::type;

	template<typename T=void>
	using Surface_d = Surface<double>;

	template<typename T=void>
	using Surface_f = Surface<float>;

	template<typename T=void>
	using Surface_ld = Surface<long double>;

	/*
	* Curves
	*/

	template<typename T>
	struct RegularCurve_ {
		using type = dg::curve::RegularCurve<param<T>>;
	};

	template<typename T>
	using RegularCurve = typename RegularCurve_<T>::type;

	template<typename T=void>
	using RegularCurve_d = RegularCurve<double>;

	template<typename T = void>
	using RegularCurve_f = RegularCurve<float>;

	template<typename T = void>
	using RegularCurve_ld = RegularCurve<long double>;

	template<typename T>
	struct UnitSpeedCurve_ {
		using type = dg::curve::UnitCurve<param<T>>;
	};

	template<typename T>
	using UnitCurve = typename UnitSpeedCurve_<T>::type;

	template<typename T=void>
	using UnitCurve_d = UnitCurve<double>;

	template<typename T = void>
	using UnitCurve_f = UnitCurve<float>;

	template<typename T = void>
	using UnitCurve_ld = UnitCurve < long double > ;

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

		template<typename T>
		struct number_helper {
			using type = boost::math::differentiation::autodiff_v1::autodiff_fvar<T, 3>;
		};

		template< typename T>
		using num = typename number_helper<T>::type;

		template<typename T>
		using nums = num<num<T>>;

	}

} //namespace dg


#endif //!TYPE_TRAITS_H
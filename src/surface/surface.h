/** @file surface.h
 *	@ingroup Surfaces
 *  @brief Smooth 2D surface class
 *
 *	A class representing a smooth two-dimensional
 *  surface parametrised by three smooth functions
 *  in two variables. It again uses special template 
 *  parameters that allow autodifferentiation to work
 *	- more details in the data members section.
 *
 *  @author Pedro Juarez Rossello (pedrojuarezrossello)
 *  @bug No known bugs.
 */

#ifndef DG_SURFACE_H
#define DG_SURFACE_H
#include "../utils/real_function.h"
#include "../vector/vector.h"
#include "../surface/plane.h"
#include "../forms/first_fundamental_form.h"
#include "../forms/second_fundamental_form.h"

namespace dg {

	namespace surf {

		/**
		* @headerfile surface.h "src/surface/surfce.h"
		* 
		* @tparam _Ret  The library has defined wrapper types 
		* around a floating point (underlying) type U
		* to allow autodifferentiation - you don't need to specify
		* this parameter. More info on the example sections and type_traits.h.
		* 
		* @tparam ...Args The library has defined two types
		* dg::param1<U> and dg::param2<U> for the first and second
		* parameter of the surface's component functions (respectively),
		* where U is a floating point (underlying) type.
		* Since the component functions have to be a template, these have to be passed
		* to the template arguments when instantiating the surface.
		*/
		template<typename _Ret, typename... Args>
		class Surface 
		{
			/**
			* @brief X component 
			*/
			dg::math::Function<_Ret, Args...> X_;

			/**
			* @brief Y component
			*/
			dg::math::Function<_Ret, Args...> Y_;

			/**
			* @brief Z component
			*/
			dg::math::Function<_Ret, Args...> Z_;

			//Derivative util function (returns a Vector)
			template<int Order1, int Order2, typename... Us>
			auto derivative_(Us... var) const
			{
				using Tuple = std::tuple<Us...>;
				//The first element of the tuple determines the type of Vector
				dg::vector::Vector< std::tuple_element<0, Tuple>::type> componentwiseDerivative( 
					X_.derivative<Order1, Order2>(var...),
					Y_.derivative<Order1, Order2>(var...),
					Z_.derivative<Order1, Order2>(var...)
				);

				return componentwiseDerivative;
			}

			//First fundamental form: E dx^2 + 2F dx dy + G dy^2

			//Computes F pointwise
			template<int Order, typename V, typename W> 
			W firstFundamentalForm__Coefficient_(V var1, V var2) const
			{
				return this->derivative_<Order, 0>(var1, var2) * this->derivative_<0, Order>(var1, var2);
			}

			//Computes E and G pointwise (Order1 = 1 and Order = 0 for E and vice versa for G)
			template<int Order1, int Order2, typename V, typename W> //E,G
			W firstFundamentalForm__Coefficient_(V var1, V var2) const
			{
				return this->derivative_<Order1, Order2>(var1, var2)* this->derivative_<Order1, Order2>(var1, var2);
			}

			//Second fundamental form: L dx^2 + 2M dx dy + N dy^2

			//Computes L, M, N pointwise depending on the values of Order1 and Order2 e.g M is Order1=Order2=1
			template<int Order1, int Order2, typename V, typename W> 
			W secondFundamentalForm__Coefficient_(V var1, V var2) const
			{
				return this->derivative_<Order1, Order2>(var1, var2)* this->normal(var1, var2);
			}

			//Pointwise Gaussian curvature 
			template<typename V>
			auto compute__GuassianCurvature_(V E, V F, V G, V L, V M, V N) const
			{
				return (L * N - M * M) / (E * G - F * F);
			}

			//Pointwise Mean curvature
			template<typename V>
			auto compute__MeanCurvature_(V E, V F, V G, V L, V M, V N) const
			{
				return (L * G - 2.0 * M * F + N * E) / (2.0 * (E * G - F * F));
			}

			//Returns E as a function
			template<typename T>
			std::function<T(T, T)> getE_() const
			{
				auto E = [this](T var1, T var2) -> T {
					return this->firstFundamentalForm__Coefficient_<1, 0, T, T>(var1, var2);
				};

				return E;
			}

			//Returns F as a function
			template<typename T>
			std::function<T(T, T)> getF_() const
			{
				auto F = [this](T var1, T var2) -> T {
					return this->firstFundamentalForm__Coefficient_<1, T, T>(var1, var2);
				};

				return F;
			}

			//Returns G as a function
			template<typename T>
			std::function<T(T, T)> getG_() const
			{
				auto G = [this](T var1, T var2) -> T {
					return this->firstFundamentalForm__Coefficient_<0, 1, T, T>(var1, var2);
				};

				return G;
			}

			//Returns L as a function
			template<typename T>
			std::function<T(T, T)> getL_() const
			{
				auto L = [this](T var1, T var2) -> T {
					return this->secondFundamentalForm__Coefficient_<2, 0, T, T>(var1, var2);
				};

				return L;
			}

			//Returns M as a function
			template<typename T>
			std::function<T(T, T)> getM_() const
			{
				auto M = [this](T var1, T var2) {
					return this->secondFundamentalForm__Coefficient_<1, 1, T, T>(var1, var2);
				};

				return M;
			}

			//Returns N as a function
			template<typename T>
			std::function<T(T, T)> getN_() const
			{
				auto N = [this](T var1, T var2) {
					return this->secondFundamentalForm__Coefficient_<0, 2, T, T>(var1, var2);
				};

				return N;
			}

			//Derivative function util
			template<int Order1, int Order2, typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			std::function<dg::vector::Vector<V>(V, V)> derivativeFunction_() const
			{
				auto derivativeLambda = [this](V var1, V var2) {
					return this->derivative_<Order1, Order2>(var1, var2);
				};

				return derivativeLambda;
			}

			//Gauss map util
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			std::function<dg::vector::Vector<V>(V, V)> gaussMap_() const
			{
				auto gaussMapLambda = [this](V var1, V var2) {
					return this->normal(var1, var2);
				};

				return gaussMapLambda;
			}

		public:

			/**
			* @brief Getter function for the X component.
			*
			* @return `const`-reference to the function.
			*/
			const auto& getX() const { return X_; }

			/**
			* @brief Getter function for the Y component.
			*
			* @return `const`-reference to the function.
			*/
			const auto& getY() const { return Y_; }

			/**
			* @brief Getter function for the Z component.
			*
			* @return `const`-reference to the function.
			*/
			const auto& getZ() const { return Z_; }

			//Constructors

			Surface() = default;

			Surface(const dg::math::Function<_Ret, Args...>& X, const dg::math::Function<_Ret, Args...>& Y, const dg::math::Function<_Ret, Args...>& Z)
				: X_(X), Y_(Y), Z_(Z) {}

			/**
			* @brief Evalutation at parameter values `...var`.
			*
			* @details It returns the result as a Vector of type V.
			*
			* @tparam ...Us Floating point types of the parameter values.
			* @param ...var Floating point values at which we evaluate the surface.
			*
			* @return Vector<T>, where T is the type of the first parameter in `...var`.
			*/
			template<typename... Us>
			auto operator()(Us ... var) 
			{
				using Tuple = std::tuple<Us...>;
				dg::vector::Vector< std::tuple_element<0, Tuple>::type> evaluation(X_(var...), Y_(var...), Z_(var...));
				
				return evaluation;
			}

			/**
			* @brief Reparametrise the curve.
			*
			* @details The function argument has to be a template function
			* that's passed template parameter `dg::param1<V>` and `dg::param2<V>`
			* for the desired underlying floating point type.
			* For instance, if we wanted to use the reparametrisation (x*y)/2
			* then we'd have to define 
			*	``template<typename T, typename U>
			*	auto repam(const T& x, const U& y) 
			*	{return 0.5*x*y;}``,
			* and then pass `repam<dg::param1_d<>, dg::param2_d<>>` to the method.
			*
			* @tparam Functor1 Function type of the first method parameter.
			* @tparam Functor2 Function type of the second method parameter.
			* @param reparametrisation1 Function representing the desired reparametrisation for the first parameter.
			* @param reparametrisation2 Function representing the desired reparametrisation for the second parameter.
			* It needs to follow the rules above.
			*
			* @return void.
			*/
			template<typename Functor1, typename Functor2>
			void reparametrise(Functor1&& reparametrisation1,Functor2&& reparametrisation2)
			{
				this->X_ = dg::math::compose<Functor1, Functor2, _Ret, Args...>(this->X_,
					std::forward<Functor1>(reparametrisation1),
					std::forward<Functor2>(reparametrisation2));
				this->Y_ = dg::math::compose<Functor1, Functor2, _Ret, Args...>(this->Y_,
					std::forward<Functor1>(reparametrisation1),
					std::forward<Functor2>(reparametrisation2));
				this->Z_ = dg::math::compose<Functor1, Functor2, _Ret, Args...>(this->Z_,
					std::forward<Functor1>(reparametrisation1), 
					std::forward<Functor2>(reparametrisation2));
			}

			/**
			* @brief Derivative function of orders Order1 and Order2.
			*
			* @details It returns an `std::function` object
			* that outputs the derivative as a Vector of type V
			* for two given values of type V.
			* Note Order1 and Order2 is a template parameter.
			*
			* @tparam Order1 Order (int) of differentiation in first component.
			* @tparam Order2 Order (int) of differentiation in second component.
			* @tparam V Floating point type of the return Vector.
			*
			* @return `std::function` containing the derivative function
			*/
			template<int Order1, int Order2, typename V>
			std::function<dg::vector::Vector<V>(V,V)> derivative() const
			{
				return derivativeFunction_<Order1,Order2,V>();
			}

			/**
			* @brief Tangent vectors at parameter values `...var`.
			* 
			* @details It computes the tangent vectors in each direction.
			* 
			* @tparam ...Us Floating point types of the parameters.
			* @param ...var Parameter values at which we compute the tangent vectors.
			* 
			* @return `std::pair` where the first is the tangent vector in the direction
			* of the first parameter and the second in the direction of the second parameter.
			*/
			template<typename... Us>
			auto tangentVectors(Us ... var) const
			{
				auto tangentVecs = std::make_pair(derivative_<1, 0>(var...), derivative_<0, 1>(var...));
				return tangentVecs;
			}

			/**
			* @brief Tangent plane at parameter values `...var`.
			* 
			* @details The tangent plane is spanned by the vectors returned from `tangentVectors`
			*
			* @tparam ...Us Floating point types of the parameters.
			* @param ...var Parameter values at which we compute the tangent plane.
			*
			* @return A Plane object.
			*/
			template<typename... Us>
			auto tangentPlane(Us ... var) 
			{
				auto tangent = tangentVectors(var...);
				Plane<Us...> plane(tangent.first, tangent.second, this->operator()(var...));
				return plane;
			}


			/**
			* @brief Unit normal vector to the surface at parameter values `...var`.
			*
			* @tparam ...Us Floating point types of the parameters.
			* @param ...var Parameter values at which we compute the normal vector.
			*
			* @return A Vector whose type is the type of the first value in the argument pack `...var`.
			*/
			template<typename... Us>
			auto normal(Us ... var) const
			{
				auto tangent = tangentVectors(var...);
				using Tuple = std::tuple<Us...>;
				dg::vector::Vector< std::tuple_element<0, Tuple>::type> normal(dg::vector::cross_product(tangent.first, tangent.second));
				normal.normalise();
				return normal;
			}

			/**
			* @brief Gaussian curvature at parameter values `...var`.
			*
			* @tparam ...Us Floating point types of the parameters.
			* @param ...var Parameter values at which we compute the Gaussian curvature.
			*
			* @return A scalar.
			*/
			template<typename... Us>
			auto gaussianCurvature(Us ... vars) const
			{
				auto L{ secondFundamentalForm__Coefficient_<2,0,Us...>(vars...) };
				auto M{ secondFundamentalForm__Coefficient_<1,1,Us...>(vars...) };
				auto N{ secondFundamentalForm__Coefficient_<0,2,Us...>(vars...) };

				auto E{ firstFundamentalForm__Coefficient_<1,0,Us...>(vars...) };
				auto F{ firstFundamentalForm__Coefficient_<1,Us...>(vars...) };
				auto G{ firstFundamentalForm__Coefficient_<0,1,Us...>(vars...) };

				return (L * N - M * M) / (E * G - F * F);
			}

			/**
			* @brief Gauss map function.
			*
			* @details It returns an `std::function` object
			* that outputs the principal normal as a Vector of type V
			* for two given values of type V.
			*
			* @tparam V Floating point type of the return Vector.
			*
			* @return `std::function` containing the Gauss map function.
			*/
			template<typename V>
			std::function<dg::vector::Vector<V>(V, V)> gaussMap() const
			{
				return gaussMap_<V>();
			}

			/**
			* @brief Mean curvature at parameter values `...var`.
			*
			* @tparam ...Us Floating point types of the parameters.
			* @param ...var Parameter values at which we compute the mean curvature.
			*
			* @return A scalar.
			*/
			template<typename... Us>
			auto meanCurvature(Us... vars) const
			{
				auto L{ secondFundamentalForm__Coefficient_<2,0,Us...>(vars...) };
				auto M{ secondFundamentalForm__Coefficient_<1,1,Us...>(vars...) };
				auto N{ secondFundamentalForm__Coefficient_<0,2,Us...>(vars...) };

				auto E{ firstFundamentalForm__Coefficient_<1,0,Us...>(vars...) };
				auto F{ firstFundamentalForm__Coefficient_<1,Us...>(vars...) };
				auto G{ firstFundamentalForm__Coefficient_<0,1,Us...>(vars...) };

				return (L * G - 2.0 * M * F + N * E) / (2.0 * (E * G - F * F));
			}

			/**
			* @brief Gaussian curvature function.
			*
			* @details It returns an `std::function` object
			* that outputs the Gaussian curvature as a scalar of type W
			* for two given values of type V.
			*
			* @tparam V Floating point type of the parameters.
			* @tparam W Floating point type of the return.
			*
			* @return `std::function` containing the Gaussian curvature function.
			*/
			template<typename V, typename W>
			std::function<W(V, V)> gaussianCurvature() const
			{
				auto gaussCurv = [this](V var1, V var2) {
					return this->gaussianCurvature(var1, var2);
				};
				return gaussCurv;
			}

			/**
			* @brief Mean curvature function.
			*
			* @details It returns an `std::function` object
			* that outputs the mean curvature as a scalar of type W
			* for two given values of type V.
			*
			* @tparam V Floating point type of the parameters.
			* @tparam W Floating point type of the return.
			*
			*
			* @return `std::function` containing the mean curvature function.
			*/
			template<typename V, typename W>
			std::function<W(V, V)> meanCurvature() const
			{
				auto meanCurv = [this](V var1, V var2) {
					return this->meanCurvature(var1, var2);
				};
				return meanCurv;
			}

			/**
			* @brief Principal curvatures at parameter values `...var`.
	
			* @tparam ...Us Floating point types of the parameters.
			* @param ...vars Parameter values at which we compute the principal curvatures.
			* 
			* @return `std::pair` where first is the maximum principal
			* curvature and second is the minimum principal curvature.
			*/
			template<typename... Us>
			auto principalCurvatures(Us... vars) const
			{
				using std::sqrt;

				auto L{ secondFundamentalForm__Coefficient_<2,0,Us...>(vars...) };
				auto M{ secondFundamentalForm__Coefficient_<1,1,Us...>(vars...) };
				auto N{ secondFundamentalForm__Coefficient_<0,2,Us...>(vars...) };

				auto E{ firstFundamentalForm__Coefficient_<1,0,Us...>(vars...) };
				auto F{ firstFundamentalForm__Coefficient_<1,Us...>(vars...) };
				auto G{ firstFundamentalForm__Coefficient_<0,1,Us...>(vars...) };

				auto K{ compute__GuassianCurvature_(E,F,G,L,M,N) };
				auto H{ compute__MeanCurvature_(E,F,G,L,M,N) };

				return std::make_pair(H + sqrt(H * H - K), H - sqrt(H * H - K));
			}

			/**
			* @brief Maximum principal curvature function.
			*
			* @details It returns an `std::function` object
			* that outputs the maximum principal curvature as a scalar of type W
			* for two given values of type V.
			*
			* @tparam V Floating point type of the parameters.
			* @tparam W Floating point type of the return.
			*
			*
			* @return `std::function` containing the maximum principal curvature function.
			*/
			template<typename V, typename W>
			std::function<W(V, V)> maxPrincipalCurvature() const
			{
				auto maxPrincipalCurv = [this](V var1, V var2) {
					return this->principalCurvatures(var1, var2).first;
				};
				return maxPrincipalCurv;
			}

			/**
			* @brief Minimum principal curvature function.
			*
			* @details It returns an `std::function` object
			* that outputs the minimum principal curvature as a scalar of type W
			* for two given values of type V.
			*
			* @tparam V Floating point type of the parameters.
			* @tparam W Floating point type of the return.
			*
			* @return `std::function` containing the minimum principal curvature function.
			*/
			template<typename V, typename W>
			std::function<W(V, V)> minPrincipalCurvature() const
			{
				auto minPrincipalCurv = [this](V var1, V var2) {
					return this->principalCurvatures(var1, var2).second;
				};
				return minPrincipalCurv;
			}

			/**
			* @brief Classify points on the surface according to curvature.
			* 
			* @details Points can be classified according to the signs of 
			* principal curvatures at that point. This method returns 
			* the classification as an enum value.
			* 
			* @tparam ...Us Floating point types of the parameters.
			* @param ...vars Parameter values at which we classify the resulting point.
			* 
			* @return Classification as an enum value Point.
			*/
			template<typename... Us>
			auto classifyPoint(Us... vars) const
			{
				auto principalCurvs = principalCurvatures(vars...);
				auto k1{ principalCurvs.first };
				auto k2{ principalCurvs.second };

				auto K{ k1 * k2 };

				if (K > 0) {
					if (!dg::math::almostEqualRelativeAndAbs(k1, k2, 1.0e-8)) return Point::elliptic;
					else return Point::umbilic;
				}
				else if (K < 0) return Point::hyperbolic;
				else {
					if ((!dg::math::isZero(k1) && dg::math::isZero(k2)) || (dg::math::isZero(k1) && !dg::math::isZero(k2))) return Point::parabolic;
					else return Point::planar;
				}
			}

			/**
			* @brief Get the first fundamental form of the surface
			* 
			* @tparam U Underlyng floating point type of the two-form.
			* 
			* @return An `std::unique_ptr` to a newly created FirstFundamentalForm object.
			*/
			template<typename U>
			std::unique_ptr<dg::form::FirstFundamentalForm<U>> firstFundamentalForm() const
			{
				return std::make_unique< dg::form::FirstFundamentalForm<U>>(getE_<U>(), getF_<U>(), getG_<U>());
			}

			/**
			* @brief Get the second fundamental form of the surface
			*
			* @tparam U Underlyng floating point type of the two-form.
			*
			* @return An `std::unique_ptr` to a newly created SecondFundamentalForm object.
			*/
			template<typename U>
			std::unique_ptr<dg::form::SecondFundamentalForm<U>> secondFundamentalForm() const
			{
				return std::make_unique< dg::form::SecondFundamentalForm<U>>(getL_<U>(), getM_<U>(), getN_<U>());
			}

		};

	} //namespace surf

} //namespace dg

#endif // !DG_SURFACE_H


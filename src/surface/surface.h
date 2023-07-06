#ifndef SURFACE_H
#define SURFACE_H
#include "../utils/real_function.h"
#include "../vector/vector.h"
#include "../surface/plane.h"
#include "../forms/first_fundamental_form.h"
#include <iostream>

namespace dg {

	namespace surf {

		template<typename _Ret, typename... Args>
		class Surface 
		{
			/*template<typename T>
			friend class TwoFormInterface;

			template<typename T>
			friend class FirstFundamentalForm;

			template<typename T>
			friend class SecondFundamentalForm;*/

			dg::math::Function<_Ret, Args...> X_;
			dg::math::Function<_Ret, Args...> Y_;
			dg::math::Function<_Ret, Args...> Z_;

			template<int Order1, int Order2, typename... Us>
			auto derivative_(Us... var) const
			{
				using Tuple = std::tuple<Us...>;
				dg::vector::Vector< std::tuple_element<0, Tuple>::type> componentwiseDerivative(
					X_.derivative<Order1, Order2>(var...),
					Y_.derivative<Order1, Order2>(var...),
					Z_.derivative<Order1, Order2>(var...)
				);

				return componentwiseDerivative;
			}


			template<int Order, typename V, typename W> //F
			W firstFundamentalForm__Coefficient_(V var1, V var2) const
			{
				return this->derivative_<Order, 0>(var1, var2) * this->derivative_<0, Order>(var1, var2);
			}

			template<int Order1, int Order2, typename V, typename W> //E,G
			W firstFundamentalForm__Coefficient_(V var1, V var2) const
			{
				return this->derivative_<Order1, Order2>(var1, var2)* this->derivative_<Order1, Order2>(var1, var2);
			}

			template<int Order1, int Order2, typename V, typename W> //F
			W secondFundamentalForm__Coefficient_(V var1, V var2) const
			{
				return this->derivative_<Order1, Order2>(var1, var2)* this->normal(var1, var2);
			}

			template<typename V>
			auto compute__GuassianCurvature_(V E, V F, V G, V L, V M, V N) const
			{
				return (L * N - M * M) / (E * G - F * F);
			}

			template<typename V>
			auto compute__MeanCurvature_(V E, V F, V G, V L, V M, V N) const
			{
				return (L * G - 2.0 * M * F + N * E) / (2.0 * (E * G - F * F));
			}

		

			template<typename ... Us>
			auto getL_() const
			{
				auto E = [this](auto... vars) {
					return this->secondFundamentalForm__Coefficient_<2, 0, Us...>(vars...);
				};

				return E;
			}

			template<typename ... Us>
			auto getM_() const
			{
				auto E = [this](auto... vars) {
					return this->secondFundamentalForm__Coefficient_<1, 1, Us...>(vars...);
				};

				return E;
			}

			template<typename ... Us>
			auto getN_() const
			{
				auto E = [this](auto... vars) {
					return this->secondFundamentalForm__Coefficient_<0, 2, Us...>(vars...);
				};

				return E;
			}


		public:

			auto const getX() const { return X_; }
			auto const getY() const { return Y_; }
			auto const getZ() const { return Z_; }

			Surface() {}

			Surface(const dg::math::Function<_Ret, Args...>& X, const dg::math::Function<_Ret, Args...>& Y, const dg::math::Function<_Ret, Args...>& Z)
				: X_(X), Y_(Y), Z_(Z) {}

			template<typename T>
			std::function<T(T, T)> getE_() const
			{
				auto E = [this](T var1, T var2) -> T {
					return this->firstFundamentalForm__Coefficient_<1, 0, T,T>(var1, var2);
				};

				return E;
			}

			template<typename T>
			std::function<T(T, T)> getF_() const
			{
				auto E = [this](T var1, T var2) -> T {
					return this->firstFundamentalForm__Coefficient_<1,T,T>(var1, var2);
				};

				return E;
			}

			template<typename T>
			std::function<T(T, T)> getG_() const
			{
				auto E = [this](T var1, T var2) {
					return this->firstFundamentalForm__Coefficient_<1, 0, T, T>(var1, var2);
				};

				return E;
			}

			//evaluation
			template<typename... Us>
			auto operator()(Us ... var) 
			{
				using Tuple = std::tuple<Us...>;
				dg::vector::Vector< std::tuple_element<0, Tuple>::type> evaluation(X_(var...), Y_(var...), Z_(var...));
				
				return evaluation;
			}

			//reparametrise
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

			//derivative function
			template<int Order1, int Order2, typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			std::function<dg::vector::Vector<V>(V,V)> derivative() const
			{
				auto derivativeLambda = [this](V var1, V var2) {
					return this->derivative_<Order1, Order2>(var1, var2);
				};

				return derivativeLambda;
			}

			//tangent vectors to a poinr
			template<typename... Us>
			auto tangentVectors(Us ... var) const
			{
				auto tangentVecs = std::make_pair(derivative_<1, 0>(var...), derivative_<0, 1>(var...));
				return tangentVecs;
			}

			//tangent space
			template<typename... Us>
			auto tangentPlane(Us ... var) 
			{
				auto tangent = tangentVectors(var...);
				Plane<Us...> plane(tangent.first, tangent.second, this->operator()(var...));
				return plane;
			}

			//standard normal
			template<typename... Us>
			auto normal(Us ... var) const
			{
				auto tangent = tangentVectors(var...);
				using Tuple = std::tuple<Us...>;
				dg::vector::Vector< std::tuple_element<0, Tuple>::type> normal(dg::vector::cross_product(tangent.first, tangent.second));
				normal.normalise();
				return normal;
			}

			//normal vector field
			template<typename V,
				std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr >
			std::function<dg::vector::Vector<V>(V, V)> gaussMap() const
			{
				auto gaussMapLambda = [this](V var1, V var2) {
					return this->normal(var1, var2);
				};

				return gaussMapLambda;
			}

			//Gaussian curvature
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

			//Mean curvature
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

			//Gaussian curvature function
			template<typename V, typename W>
			std::function<W(V, V)> gaussianCurvature() const
			{
				auto gaussCurv = [this](V var1, V var2) {
					return this->gaussianCurvature(var1, var2);
				};
				return gaussCurv;
			}

			//Mean curvature function
			template<typename V, typename W>
			std::function<W(V, V)> meanCurvature() const
			{
				auto meanCurv = [this](V var1, V var2) {
					return this->meanCurvature(var1, var2);
				};
				return meanCurv;
			}

			//Principal curvatures
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

			//Maximum principal curvature function
			template<typename V, typename W>
			std::function<W(V, V)> maxPrincipalCurvature() const
			{
				auto maxPrincipalCurv = [this](V var1, V var2) {
					return this->principalCurvatures(var1, var2).first;
				};
				return maxPrincipalCurv;
			}

			//Minimum principal curvature function
			template<typename V, typename W>
			std::function<W(V, V)> minPrincipalCurvature() const
			{
				auto minPrincipalCurv = [this](V var1, V var2) {
					return this->principalCurvatures(var1, var2).second;
				};
				return minPrincipalCurv;
			}

			//Classify a point according to curvature
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

			template<typename U>
			dg::form::FirstFundamentalForm<U> getFFF() const
			{
				dg::form::FirstFundamentalForm<U> fff(getE_<U>(), getF_<U>(), getG_<U>());
				return fff;
			}

		};

	} //namespace surf

} //namespace dg

#endif // !SURFACE_H


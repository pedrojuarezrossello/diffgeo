#ifndef MATHS_H
#define MATHS_H
#include <cfloat>
#include <cmath>
#include <functional>
#include "../utils/real_function.h"

namespace dg {

    namespace math {

        constexpr auto PI = 3.14159265358979323846;

        template<typename T, typename U,
            std::enable_if_t<std::is_floating_point_v<T>&& std::is_floating_point_v<U>, std::nullptr_t> = nullptr>
        bool almostEqualRelativeAndAbs(T A, U B,
            double maxDiff, float maxRelDiff = FLT_EPSILON)
        {
            // Check if the numbers are really close -- needed
            // when comparing numbers near zero.
            auto diff = fabs(A - B);
            if (diff <= maxDiff)
                return true;

            A = fabs(A);
            B = fabs(B);

            auto largest = (B > A) ? B : A;

            if (diff <= largest * maxRelDiff)
                return true;
            return false;
        }

        template<typename T,
            std::enable_if_t<std::is_floating_point_v<T>, std::nullptr_t> = nullptr>
        bool isZero(T a)
        {
            return almostEqualRelativeAndAbs(a, 0.0, 1.0e-14);
        }

        template<typename T,
            std::enable_if_t<std::is_floating_point_v<T>, std::nullptr_t> = nullptr>
        auto signOf(T num)
        {
            return (num >= 0) ? "+" : "";
        }

     
        template <typename Functor, typename _Ret, typename Arg> 
        auto compose(const dg::math::Function<_Ret,Arg>& f, const Functor & g) 
        {
            dg::math::Function<_Ret, Arg> compose_([=](auto x) {
                return f(g(x));
                });

            return compose_;
        }

       template <typename Functor, typename _Ret, typename Arg>
        auto compose(const dg::math::Function<_Ret,Arg>& f, Functor&& g)
        {
            dg::math::Function<_Ret, Arg> compose_([=](auto x) {
                return f(g(x));
                });

            return compose_;
        }

        template <typename Functor1,typename Functor2, typename _Ret, typename ... Args>
        auto compose(const dg::math::Function<_Ret, Args...>& f, Functor1&& g, Functor2&& h)
        {
            dg::math::Function<_Ret, Args...> compose_([=](auto... x) {
                return f(g(x...), h(x...));
                    });
        
            return compose_;
        }

    } //namespace math

    namespace surf
    {
        enum class Point
        {
            elliptic,
            hyperbolic,
            parabolic,
            planar,
            umbilic
        };

        constexpr const char* typeToString(Point point) 
        {
            switch (point)
            {
            case Point::elliptic: return "Elliptic";
            case Point::parabolic: return "Parabolic";
            case Point::hyperbolic: return "Hyperbolic";
            case Point::planar: return "Planar";
            case Point::umbilic: return "Umbilic";
            default: return "Non-computable";
            }
        }

    } //namespace surf

} //namespace dg

#endif //!MATHS_H
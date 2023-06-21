#ifndef CURVE_H
#define CURVE_H
#include <type_traits>
#include "../utils/real_function.h"

template<typename T,
	std::enable_if_t<std::is_floating_point_v<T>, std::nullptr_t> = nullptr>
class Curve 
{
	Function __X;
	Function __Y;
	Function __Z;

public:
	//double operator()(T var);
};

#endif //!CURVE_H

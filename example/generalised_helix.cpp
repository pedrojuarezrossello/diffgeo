#include <iostream>
#include "src/curve/curve.h"

/*PROBLEM STATEMENT:

Given the curve a(u) = (2 cos u + 2 sin u, 2 sin u − 2 cos u, u),
(a)find the length between p = a(0) and u = a(2*pi)
(b)Show that a is a generalised helix, i.e. the tangents to a make a constant angle with some
fixed direction
(c)Compute curvature and torsion 
*/

template<typename T>
auto helix1(T const& x) {
	using std::cos;
	return 2.0 * (cos(x) + sin(x));
}

template<typename T>
auto helix2(T const& x) {
	using std::sin;
	return 2.0 * (sin(x)-cos(x));
}

template<typename T>
auto helix3(T const& x) {
	return x;
}

//SET UP
double pi = dg::math::PI;
dg::CurveComponent_d<> X_2(helix1<dg::param_d<>>);
dg::CurveComponent_d<> Y_2(helix2<dg::param_d<>>);
dg::CurveComponent_d<> Z_2(helix3<dg::param_d<>>);
dg::RegularCurve_d<> helix(X_2, Y_2, Z_2);

int main()
{
	//a) Find length
	std::cout << helix.length(0.0, 2 * pi) << std::endl;
	
	//Result = 18.8496 ~ 6*pi

	//b) The tangents to a make a constant angle with some fixed direction

	/*We can't exactly prove it, but we can see that three random point 
	have all the same angle with the vector (0,0,1), as required.*/
	auto tangent1(helix.unitTangent(0.0));
	auto tangent2(helix.unitTangent(2.8));
	auto tangent3(helix.unitTangent(1.9 * pi));

	dg::vector::Vector<double> axis(0.0, 0.0, 1.0);

	//4.37255
	std::cout << angle(tangent1, axis) << " " << angle(tangent2, axis) << " " << angle(tangent2, axis) << std::endl; 

	//c) Compute curvature and torsion 
	auto curvatureFunction{ [](auto var) {
		return helix.curvature(var);
		}
	};

	auto torsionFunction{ [](auto var) {
		return helix.torsion(var);
		}
	};

	//Example of usage
	std::cout << "Curvature at u=1: " << curvatureFunction(1.0) << std::endl; // 0.31427
	std::cout << "Torsion at u=1: " << torsionFunction(1.0); //-0.111111 

	return 0;
}
#include <iostream>
#include "src/surface/surface.h"
#include "src/forms/first_fundamental_form.h"
#include <boost/math/quadrature/gauss.hpp>

/**
* PROBLEM STATEMENT:
* 
* Compute the arc length of a curve u=t, v=t for 0<=t<=1 on
* a hyperbolic paraboloid r(u,v)=(u,v,uv) where 0<=u,v<=1.
* 
* The exact solution is sqrt(3/2)+0.5*ln(sqrt(2)+sqrt(3)) 
* which is approximately 1.79785.
*/

//SET UP

template<typename T>
auto xCurve_(const T& t) { return t; }

template<typename T>
auto yCurve_(const T& t) { return t; }

template<typename U,typename V>
auto xParaboloid_(const U& u, const V& v) { return u; }

template<typename U, typename V>
auto yParaboloid_(const U&, const V& v) { return v; }

template<typename U, typename V>
auto zParaboloid_(const U& u, const V& v) { return u * v; }

//Library's own types compatible witg autodifferentiation
dg::CurveComponent_d<> xCurve(xCurve_<dg::param_d<>>);
dg::CurveComponent_d<> yCurve(yCurve_<dg::param_d<>>);

dg::SurfaceComponent_d<> xParaboloid(xParaboloid_<dg::param1_d<>, dg::param2_d<>>);
dg::SurfaceComponent_d<> yParaboloid(yParaboloid_<dg::param1_d<>, dg::param2_d<>>);
dg::SurfaceComponent_d<> zParaboloid(zParaboloid_<dg::param1_d<>, dg::param2_d<>>);

int main()
{
	//Define the surface
	dg::Surface_d<> paraboloid(xParaboloid, yParaboloid, zParaboloid);

	//Compute the paraboloid's first fundamental form
	auto firstFundForm{ paraboloid.firstFundamentalForm<double>() };
	
	//Get the line element along the curve components
	auto paraboloidLineElement{ firstFundForm->lineElement(xCurve,yCurve) };

	//Integrate line element 
	std::cout << boost::math::quadrature::gauss<double, 7>::integrate(paraboloidLineElement, 0.0, 1.0);

	return 0;
}
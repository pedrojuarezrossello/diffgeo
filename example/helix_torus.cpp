#include <iostream>
#include "src/curve/curve.h"
#include "src/utils/type_traits.h"
#include "src/utils/real_function.h"
#include "src/curve/unit_curve.h"
#include "src/surface/surface.h"
#include "src/forms/first_fundamental_form.h"
#include "src/forms/vector_field.h"
#include <tuple>
#include <chrono>

template<typename T>
auto helix1(T const& x) {
	using std::cos;
	return 2.0 * cos(x);
}

template<typename T>
auto helix2(T const& x) {
	using std::sin;
	return 2.0 * sin(x);
}

template<typename T>
auto helix3(T const& x) {
	return std::sqrt(8) * x;
}

template<typename T>
auto unitReparametrisation(T const& x) {
	return  x/std::sqrt(12);
}

template<typename T, typename U>
auto torus1(T const& x, U const& y) {
	using std::cos;
	return (2 + cos(y)) * cos(x);
}

template<typename T, typename U>
auto torus2(T const& x, U const& y) {
	using std::cos;
	using std::sin;
	return (2 + cos(y)) * sin(x);
}

template<typename T, typename U>
auto torus3(T const& x, U const& y) {
	using std::sin;
	return sin(y);
}

int main() {

	using namespace std::chrono;

	time_point<high_resolution_clock> start_point, end_point; 

	start_point = high_resolution_clock::now();



		//CURVE -----------------------------

			//SET UP
		auto pi = dg::math::PI;
		dg::CurveComponent_d<> X_2(helix1<dg::param_d<>>);
		dg::CurveComponent_d<> Y_2(helix2<dg::param_d<>>);
		dg::CurveComponent_d<> Z_2(helix3<dg::param_d<>>);
		dg::RegularCurve_d<> helix(X_2, Y_2, Z_2);

		//CALCULATIONS

		helix(pi);
		auto derHelix = helix.derivative<1, double>();
		derHelix(pi);
		helix.curvature(pi);
		helix.torsion(pi);
		helix.length(0.0, 2.0 * pi);
		helix.unitTangent(pi);
		auto unitHelix = helix.unitSpeedParametrisation(unitReparametrisation<dg::param_d<>>);
		unitHelix->curvature(pi);
		unitHelix->torsion(pi);
		unitHelix->length(0.0, 2.0 * pi * std::sqrt(12));
		unitHelix->unitTangent(pi);
		unitHelix->principalNormal(pi);
		unitHelix->binormal(pi);
		unitHelix->osculatingPlane(pi);
		unitHelix->rectifyingPlane(pi);
		unitHelix->normalPlane(pi);



		//SURFACE ---------------------------

		//SET UP
		dg::SurfaceComponent_d<> X(torus1<dg::param1_d<>, dg::param2_d<>>);
		dg::SurfaceComponent_d<> Y(torus2<dg::param1_d<>, dg::param2_d<>>);
		dg::SurfaceComponent_d<> Z(torus3<dg::param1_d<>, dg::param2_d<>>);
		dg::Surface_d<> torus(X, Y, Z);

		torus(pi, 2 * pi);
		auto tangVec = torus.tangentVectors(pi, 2 * pi);
		torus.tangentPlane(pi, 2 * pi);
		torus.normal(pi, 2 * pi);
		auto gaussMap = torus.gaussMap<double>();
		gaussMap(pi, 2 * pi);
		torus.gaussianCurvature(pi, 2 * pi);
		auto gaussCurv = torus.gaussianCurvature<double, double>();
		gaussCurv(pi, 2 * pi);
		torus.meanCurvature(2 * pi, pi / 2);
		auto meanCurv = torus.meanCurvature<double, double>();
		meanCurv(2 * pi, pi / 2);
		torus.principalCurvatures(pi, 2 * pi);
		auto minPrincipal = torus.minPrincipalCurvature<double, double>();
		minPrincipal(pi, 2 * pi);
		dg::surf::typeToString(torus.classifyPoint(pi, 2 * pi));
		torus.reparametrise(torus1<dg::param1_d<>, dg::param2_d<>>, torus2<dg::param1_d<>, dg::param2_d<>>);



		//FORMS --------------------------------

		//SET UP
		dg::form::VectorField<double> vectorField1(torus1<double, double>, torus2<double, double >);
		dg::form::VectorField<double> vectorField2(torus2<double, double>, torus3<double, double >);
		auto firstFF = torus.firstFundamentalForm<double>();
		auto secondFF = torus.secondFundamentalForm<double>();

		//CALCULATIONS
		auto areaEl = firstFF->areaElement();
		areaEl(pi, pi);
		auto lineEl = firstFF->lineElement<dg::param_d<>>(X_2, X_2);
		lineEl(pi);
		vectorField1.compute_U(pi, pi);
		secondFF->innerProduct<double>(vectorField1, vectorField2, pi, pi * 0.5);
	

	end_point = high_resolution_clock::now(); 

	auto start = time_point_cast<microseconds>(start_point).time_since_epoch().count();

	auto end = time_point_cast<microseconds>(end_point).time_since_epoch().count();

	std::cout << std::endl << "Time taken = " << (end - start) << " microseconds" << std::endl;

	return 0;
}
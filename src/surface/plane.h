#ifndef PLANE_H
#define PLANE_H
#include <cstddef>
#include "../vector/vector.h"
#include "../utils/operations.h"
#include "../curve/line.h"
#include "../utils/maths.h"

template<typename T, typename U,
         std::enable_if_t<std::is_floating_point_v<T>&& std::is_floating_point_v<U>, std::nullptr_t> = nullptr>
class Plane
{
    Vector<T> normal_;
    Vector<U> point_;

public:

    //constructors & equality

    Plane() = default;

    Plane(const Vector<T>& normal, const Vector<U>& point)
        : normal_(normal), point_(point)
    {
        normal_.normalise();
    }

    Plane(const Vector<T>& vec1, const Vector<T>& vec2, const Vector<U>& point)
	    : normal_(cross_product(vec1, vec2)), point_(point)
    {
        normal_.normalise();
    }

    template<typename V,
        std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr>
    Plane(const Line<T,V>& line, const Vector<U>& point)
	    : normal_(cross_product(line.getDirection(),line.at(0.0)-point)), point_(point)
    {
        normal_.normalise();
    }

    template<typename V, typename W,
        std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
    Plane(const Line<T,V>& line1, const Line<U,W>& line2)
	    : normal_(Vector<T>()), point_(Vector<U>())
    {
        if (line1.intersects(line2)) 
        {
            normal_ = cross_product(line1.getDirection(), line2.getDirection);  normal_.normalise();
            point_ = line1.at(0.0);
        }
        else if (areParallel(line1,line2))
        {
            normal_ = cross_product(line1.getDirection(), line1.at(0.0) - line2.at(0.0)); normal_.normalise();
            point_ = line1.at(0.0);
        }
        else { throw "Lines must intersect or be parallel"; }
    }

    template <typename CallableObject, typename... Args>
    Plane(const VectorExpression<CallableObject,Args...>& vector_expression, const Vector<U>& point)
	    : normal_(vector_expression), point_(point)
    {
        normal_.normalise();
    }

    template <typename CallableObject, typename... Args>
    Plane(const Vector<T>& normal, const VectorExpression<CallableObject, Args...>& point)
        : normal_(normal), point_(point)
    {
        normal_.normalise();
    }

    template <typename CallableObject1, typename... Args1, typename CallableObject2, typename... Args2>
    Plane(const VectorExpression<CallableObject1, Args1...>& vector_expression, const VectorExpression<CallableObject2, Args2...>& point)
	    : normal_(vector_expression), point_(point)
    {
        normal_.normalise();
    }

    bool operator==(const Plane& other) const
    {
        return normal_ == other.normal_ && point_ == other.point_;
    }

    //getters
    const Vector<T>& getNormal() const
    {
        return normal_;
    }

    Vector<T>& getDirection()
    {
        return const_cast<int&>(const_cast<const Plane*>(this)->getNormal());
    }

    //print util
    friend ostream& operator<<(ostream& os, Plane const& plane)
    {
    	auto sign{""};

        const T scalar{ plane.normal_ * plane.point_ };
        if (!isZero(plane.normal_[0]))
        {
            os << plane.normal_[0] << "x";
        }
        if (!isZero(plane.normal_[1]))
        {
            sign = signOf(plane.normal_[1]);
            os << sign << plane.normal_[1] << "y";
        }
        if (!isZero(plane.normal_[2]))
        {
            sign = signOf(plane.normal_[2]);
            os << sign << plane.normal_[2] << "z";
        }
        sign = signOf(scalar);
        os << sign << scalar << "=0";
        return os;
    }

    //is a point on the plane?
    template<typename V,
        std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr>
    bool isOnPlane(const Vector<V>& point)
    {
        return isZero(normal_ * point - point_);
    }

    //is a line on the plane?
    template<typename V, typename W,
        std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
    bool isOnPlane(const Line<V,W>& line)
    {
        return isOnPlane(line.at(0.0)) && isOnPlane(line.at(1.0));
    }

    //projections
    template<typename V,
        std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr>
    auto projection(const Vector<V>& point)
    {
        Vector<V> projOnNormal = Vector((point * normal_) * normal_);
        Vector<V> projOnPlane = point - projOnNormal;
        return projOnPlane;
    }

    template<typename V, typename W,
        std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
    auto projection(const Line<V,W>& line)
    {
        auto vec = Vector(this->projection(line.getDirection()));
        auto point = this->intersectsAt(line);
        auto projLine = Line(point, vec);
        return projLine;
    }

    //distance
    template<typename V,
        std::enable_if_t<std::is_floating_point_v<V>, std::nullptr_t> = nullptr>
    double distanceFrom(const Vector<V>& point)
    {
        return point * normal_;
    }

    //line & plane

    //parallel
    template<typename V, typename W,
        std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
    friend bool areParallel(const Plane& plane1, const Plane<V,W>& plane2)
    {
        return plane1.normal_ == plane2.normal_;
    }

    template<typename V, typename W,
        std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
    friend bool areParallel( Plane& plane1, const Plane<V, W>& plane2)
    {
        return plane1.normal_ == plane2.normal_;
    }

    template<typename V, typename W,
        std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
    friend bool areParallel(const Plane& plane, const Line<V, W>& line)
    {
        return plane.normal_.isPerpendicular(line.getDirection());
    }

    template<typename V, typename W,
        std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
    friend bool areParallel(Plane& plane, const Line<V, W>& line)
    {
        return plane.normal_.isPerpendicular(line.getDirection());
    }

    //intersection
    template<typename V, typename W,
        std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
    bool intersects(const Line<V, W>& line)
    {
        return !areParallel(*this, line);
    }

    template<typename V, typename W,
        std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
    bool intersects(const Plane<V, W>& plane)
    {
        return !areParallel(*this, plane);
    }

    template<typename V, typename W,
        std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
    auto intersectsAt(const Line<V, W>& line)
    {
        if (!this->intersects(line)) throw std::exception("Lines do not intersect!");

        const auto parameter{ ((point_ - line.at(0.0)) * normal_) / (line.getDirection() * normal_) };
        auto pointIntersection = Vector(line.at(parameter));
        return pointIntersection;
      
    }

    template<typename V, typename W,
        std::enable_if_t<std::is_floating_point_v<V>&& std::is_floating_point_v<W>, std::nullptr_t> = nullptr>
    auto intersectsAlong(const Plane<V, W>& plane)
    {
        if (!this->intersects(plane)) throw std::exception("Planes don't intersect!");

        const auto planeConst1{ normal_ * point_ };
        const auto planeConst2{ plane.normal_ * plane.point_ };
        const auto dottedNormals{ normal_ * plane.normal_ };

        const auto coeff1{ (planeConst1 - planeConst2 * dottedNormals) / (1-dottedNormals*dottedNormals)};
        const auto coeff2{ (planeConst2 - planeConst1 * dottedNormals) / (1 - dottedNormals * dottedNormals) };

        auto lineIntersection( Line<T,U>(coeff1 * normal_ + coeff2 * plane.normal_,cross_product(normal_,plane.normal_)) );
        return lineIntersection;

    }


};


#endif //!PLANE_H
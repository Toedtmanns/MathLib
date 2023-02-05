#pragma once
#include "mlPrimitives.hpp"
#include "mlMatrices.hpp"
#include "mlVectors.hpp"
#include <string.h>

namespace MathLib
{
	// Lines

	struct EXPORT Intersect
	{
		Float2 pos;
		bool isIntersecting, isCollinear;
		Line2D line1, line2;

		Intersect();
		Intersect(const Float2& pos, const bool intersecting, const bool collinear);
	};

	EXPORT float GetSlope(const Line2D& line);
	EXPORT float GetYAxisSection(const Line2D& line);
	EXPORT bool IsParallel(const Line2D& l1, const Line2D& l2);
	EXPORT Float2 GetIntersectPos(const Line2D& l1, const Line2D& l2);
	EXPORT Intersect GetIntersect(const Line2D& l1, const Line2D& l2);
	EXPORT Intersect GetIntersect(const Line2D& l, const Float2& p);

	EXPORT void PrintProperties(const Line2D& l);
	EXPORT void PrintProperties(const Intersect& i);

	EXPORT float Lerp(const float start, const float end, const float t);
	EXPORT Float2 Lerp(const Float2& start, const Float2& end, const float t);
	EXPORT Float2 Lerp(const Line2D& line, const float t);
	EXPORT Float3 Lerp(const Float3& start, const Float3& end, const float t);
	EXPORT Float3 Lerp(const Line3D& line, const float t);

	// Quaternions

	class EXPORT Quaternion
	{
	public:
		float real, i, j, k;

		Quaternion();
		Quaternion(const float real, const float i, const float j, const float k);
		Quaternion(const Float3& point);

		Quaternion operator+(const Quaternion& other) const;
		Quaternion operator*(const Quaternion& other) const;
		void operator+=(const Quaternion& other);
		void operator*=(const Quaternion& other);
		Quaternion RotateQuaternion(const Quaternion& quat) const;
		Float3 RotatePoint(const Float3& point) const;

		Quaternion GetInverse() const;
		Float3 GetPoint() const;
	};

	EXPORT Quaternion QuaternionRotation(const float angle, const float iAxis, const float jAxis, const float kAxis);
	EXPORT Quaternion QuaternionRotation(const float angle, const Float3& axis);
	EXPORT void PrintProperties(const Quaternion& quat);

	// Quaternion operations

	EXPORT Mat4 MatrixRotate(const Mat4& mat, const Quaternion& quat);

	// Geometry maths

	template<typename T>
	class BaseIterator
	{
	protected:
		using ValType = typename T::ValueType;
		ValType* m_Ptr;
	public:
		BaseIterator(ValType* ptr)
			: m_Ptr(ptr)
		{

		}

		BaseIterator& operator++()
		{
			m_Ptr++;
			return *this;
		}
		BaseIterator operator++(int)
		{
			BaseIterator tmp = *this;
			++(*this);
			return tmp;
		}
		BaseIterator& operator--()
		{
			m_Ptr--;
			return *this;
		}
		BaseIterator operator--(int)
		{
			BaseIterator tmp = *this;
			--(*this);
			return tmp;
		}
		ValType& operator[](const size_t index)
		{
			return *(m_Ptr + index);
		}
		ValType* operator->()
		{
			return m_Ptr;
		}
		ValType& operator*()
		{
			return *m_Ptr;
		}
		bool operator==(const BaseIterator& other) const
		{
			return m_Ptr == other.m_Ptr;
		}
		bool operator!=(const BaseIterator& other) const
		{
			return !(m_Ptr == other.m_Ptr);
		}
	};

	class EXPORT GeoCollision
	{
	protected:
		size_t m_IntersectArrLength;
		size_t m_IntersectCount;
		Intersect* m_IntersectArray;

	public:
		using Iterator = BaseIterator<GeoCollision>;
		using ValueType = Intersect;
		bool isColliding;

		GeoCollision();
		GeoCollision(GeoCollision&& other) noexcept;
		GeoCollision(const GeoCollision& other);
		GeoCollision(const Intersect* intersectArr, const size_t intersectArrLength);

		GeoCollision& AddIntersect(const Intersect& intersect);
		size_t GetIntersectCount() const;

		void operator=(GeoCollision&& other) noexcept;
		void operator=(const GeoCollision& other);
		Intersect& operator[](const size_t index);

		Iterator begin();
		Iterator end();

		~GeoCollision();
	};

	class Circle2D;

	class EXPORT Polygon2D
	{
	public:
		using Iterator = BaseIterator<Polygon2D>;
		using ValueType = Float2;
		size_t m_NumCorners;
		Float2* m_CornerArr;

		Polygon2D();
		Polygon2D(Polygon2D&& other) noexcept;
		Polygon2D(const Polygon2D& other);
		Polygon2D(const size_t corners);
		Polygon2D(const Float2* const pointArr, const size_t corners);

		Polygon2D& Translate(const Float2& translation);
		Polygon2D& Translate(const float translateX, const float translateY);
		Polygon2D& Rotate(const float angle);
		Polygon2D& Scale(const float scale);
		Polygon2D& Scale(const float scaleX, const float scaleY);

		Float2 getCenter() const;
		Float2 SupportFunction(const Vector2D& direction) const;
		bool CollidesWith(const Polygon2D& other) const;
		bool CollidesWith(const Circle2D& other) const;
		bool Contains(const Float2& point) const;

		void operator=(Polygon2D&& other) noexcept;
		void operator=(const Polygon2D& other);

		Iterator begin();
		Iterator end();

		~Polygon2D();
	};

	class EXPORT Circle2D
	{
	public:
		Float2 position;
		float radius;

		Circle2D();
		Circle2D(const Float2& position, const float radius = 0.5);

		Circle2D& Translate(const Float2& translation);
		Circle2D& Translate(const float translateX, const float translateY);

		Float2 SupportFunction(const Vector2D& direction) const;
		bool CollidesWith(const Polygon2D& other) const;
		bool CollidesWith(const Circle2D& other) const;
	};

	class EXPORT Triangle2D : public Polygon2D
	{
	public:
		Triangle2D();
		Triangle2D(Triangle2D&& other) noexcept;
		Triangle2D(const Triangle2D& other);
		Triangle2D(const Float2& p1, const Float2& p2, const Float2& p3);
		Triangle2D(const Float2* const pointArr);

		void operator=(Triangle2D&& other) noexcept;
		void operator=(const Triangle2D& other);
	};

	class EXPORT Rectangle2D : public Polygon2D
	{
	public:
		Rectangle2D();
		Rectangle2D(Rectangle2D&& other) noexcept;
		Rectangle2D(const Rectangle2D& other);
		Rectangle2D(const Float2& p1, const Float2& p2, const Float2& p3, const Float2& p4);
		Rectangle2D(const Float2* const pointArr);
		Rectangle2D(const Float2& position, const float rotation, const Float2& size);

		void operator=(Rectangle2D&& other) noexcept;
		void operator=(const Rectangle2D& other);
	};

	class EXPORT Quad
	{
	protected:
		Float2 m_P1, m_P2, m_P3, m_P4;

		void CalcP4();

	public:
		Quad();
		Quad(const Float2& p1, const Float2& p2, const Float2& p3);

		Quad& SetQuad(const Float2& p1, const Float2& p2, const Float2& p3);
		Quad& SetP1(const Float2& point);
		Quad& SetP2(const Float2& point);
		Quad& SetP3(const Float2& point);

		Float2 GetP1() const;
		Float2 GetP2() const;
		Float2 GetP3() const;
		Float2 GetP4() const;
		Float2* GetPoints() const;
		bool Contains(const Float2& point) const;
		bool Contains(const float xPos, const float yPos) const;
	};

	class EXPORT Plane
	{
	/// <summary>
	/// A class modelling a 3-dimensional plane defined by 3 points.
	/// The 3 points are expected to be in clockwise sequence.
	/// </summary>
	protected:
		Float3 m_P1, m_P2, m_P3;
		Vector3D m_Normal;

		void CalcNormal();

	public:
		Plane();
		Plane(const Float3& p1, const Float3& p2, const Float3& p3);
		Plane(const Mat3& matrix);

		Plane& SetPlane(const Float3& p1, const Float3& p2, const Float3& p3);
		Plane& SetPlane(const Mat3& matrix);
		Plane& SetP1(const Float3& point);
		Plane& SetP2(const Float3& point);
		Plane& SetP3(const Float3& point);
		Plane& Transform(const Mat3& mat);
		Plane& Translate(const Float3& translation);

		bool IsIntersecting(const Line3D& line) const;
		bool GetIntersection(const Line3D& line, Float3* intersection) const;
		//bool IsIntersecting(const Plane& plane) const;
		//Line3D GetIntersection(const Plane& plane, bool* isIntersecting) const;

		Vector3D GetNormal() const;
		Float3 GetP1() const;
		Float3 GetP2() const;
		Float3 GetP3() const;
		Float3 GetP4() const;
		Mat3 GetAsMat() const;
	};

	EXPORT bool Contains(const Rectangle2D& rect, const float rotation, const Float2& point);
	EXPORT float* ProjectTo1D(const Polygon2D& polygon, const Vector2D& viewDir);
	EXPORT float* ProjectTo1D(const Float2* pointArr, const size_t pointCount, const Vector2D& viewDir);

	// Utility
	template <typename T>
	T MinFromArray(const T* arr, const size_t length)
	{
		T Min = arr[0];
		for (size_t i = 1; i < length; i++)
		{
			if (arr[i] < Min)
				Min = arr[i];
		}
		return Min;
	}
	template <typename T>
	T MaxFromArray(const T* arr, const size_t length)
	{
		T Max = arr[0];
		for (size_t i = 1; i < length; i++)
		{
			if (arr[i] > Max)
				Max = arr[i];
		}
		return Max;
	}
	EXPORT Line2D Vector2Line(Vector2D vector, Float2 pos);
	EXPORT int RandInt(int Min, int Max);
	inline double RoundFTo(const double num, const int decimal)
	{
		return round(Pow(10, decimal) * num) / Pow(10, decimal);
	}
	inline float RoundFTo(const float num, const int decimal)
	{
		return roundf(Pow(10, decimal) * num) / Pow(10, decimal);
	}

	template <typename T>
	void SetArraySize(T** arr, const size_t& curLen, const size_t& newLen)
	{
		T* newArr = new T[newLen];
		memcpy(newArr, *arr, Min(curLen, newLen) * sizeof(**arr));
		delete[] *arr;
		*arr = newArr;
	}

	EXPORT void PrintProperties(const Float2& p);
	EXPORT void PrintProperties(const Float3& p);
}
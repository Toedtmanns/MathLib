#pragma once
#include "mlPrimitives.hpp"
#include "mlMatrices.hpp"
#include <string.h>

namespace MathLib
{
	EXPORT void PrintProperties(const Float2& p);
	EXPORT void PrintProperties(const Float3& p);

	// Lines

	struct EXPORT Line2D
	{
		Float2 p1;
		Float2 p2;
		double normal;

		Line2D();
		Line2D(const float x1, const float y1, const float x2, const float y2);
		Line2D(const Float2& p1, const Float2& p2);

		bool operator==(const Line2D& other);
		bool operator!=(const Line2D& other);
	};

	struct EXPORT Intersect
	{
		Float2 pos;
		bool isIntersecting, isCollinear;
		Line2D line1, line2;

		Intersect();
		Intersect(const Float2& pos, const bool intersecting, const bool collinear);
	};

	EXPORT float GetLength(const Line2D& line);
	EXPORT float GetDistance(const Float2& p1, const Float2& p2);
	EXPORT float GetDistance(const Float3& p1, const Float3& p2);

	EXPORT float GetSlope(const Line2D& line);
	EXPORT float GetYAxisSection(const Line2D& line);
	EXPORT bool IsParallel(const Line2D& l1, const Line2D& l2);
	EXPORT Float2 GetIntersectPos(const Line2D& l1, const Line2D& l2);
	EXPORT Intersect GetIntersect(const Line2D& l1, const Line2D& l2);
	EXPORT Intersect GetIntersect(const Line2D& l, const Float2& p);

	EXPORT void PrintProperties(const Line2D& l);
	EXPORT void PrintProperties(const Intersect& i);

	float Lerp(const float start, const float end, const float t);
	Float2 Lerp(const Float2& start, const Float2& end, const float t);
	Float2 Lerp(const Line2D& line, const float t);

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

	// Vector math

	class EXPORT Vector2D : public Float2
	{
	public:
		Vector2D();
		Vector2D(const Float2& dir);
		explicit Vector2D(const float x, const float y);
		Vector2D(const Line2D& line);
		Vector2D(const Float2& p1, const Float2& p2);

		const Vector2D& Transform(const Mat2& transformMat);
		const Vector2D& Scale(const float scale);
		const Vector2D& Scale(const float scaleX, const float scaleY);
		const Vector2D& SetScale(const float scale);
		const Vector2D& Rotate(float angle);
		const Vector2D& Rot90R();
		const Vector2D& Rot90L();
		const Vector2D& Normalize();

		Vector2D TripleProduct(const Vector2D& other) const;

		Vector2D operator+(const Vector2D& other) const;
		Vector2D operator-(const Vector2D& other) const;
		float operator*(const Vector2D& other) const;
		Vector2D operator*(const float number) const;
		Vector2D operator*(const Mat2& matrix) const;
		void operator*=(const float number);
		Vector2D operator-() const;

		float GetAngle() const;
		float GetLen() const;
		Mat2x1 GetRowVector() const;
		Mat1x2 GetColVector() const;

		Float2 TransformPoint(Float2 point) const;
	};

	class EXPORT Vector3D : public Float3
	{
	public:
		Vector3D();
		Vector3D(const Float3& dir);
		explicit Vector3D(const float x, const float y, const float z);
		Vector3D(const Vector2D& vec, const float z = 0);
		Vector3D(const Float3& p1, const Float3& p2);

		const Vector3D& Transform(const Mat3& mat);
		const Vector3D& Scale(const float& scale);
		const Vector3D& Scale(const float& scaleX, const float& scaleY, const float& scaleZ);
		const Vector3D& SetScale(const float& len);
		const Vector3D& Rotate(const float& angle, const unsigned int& axis);
		const Vector3D& Normalize();

		Vector3D operator+(const Vector3D& other) const;
		Vector3D operator-(const Vector3D& other) const;
		float operator*(const Vector3D& other) const;
		Vector3D operator*(const float number) const;
		Vector3D operator*(const Mat3& matrix) const;
		void operator*=(const float number);
		Vector3D operator-() const;

		float GetAngle(const int axis) const;
		float GetLen() const;
		Mat3x1 GetRowVector() const;
		Mat1x3 GetColVector() const;

		Float3 TransformPoint(Float3 point);
	};

	EXPORT Vector2D operator*(const Mat2 & matrix, const Vector2D& vector);
	EXPORT Vector3D operator*(const Mat3& matrix, const Vector3D& vector);

	EXPORT Vector3D GetRelativeVec(const Vector3D& vec1, const Vector3D& vec2);
	EXPORT float VectorGetAngleDifference(const Vector2D& v1, const Vector2D& v2);
	EXPORT float VectorDotProduct(const Vector2D& v1, const Vector2D& v2);
	EXPORT Vector3D VectorCrossProduct(const Vector2D& v1, const Vector2D& v2);
	EXPORT Vector3D VectorCrossProduct(const Vector3D& v1, const Vector3D& v2);
	EXPORT Vector2D VectorTripleProduct(const Vector2D& v1, const Vector2D& v2);
	EXPORT Vector3D VectorTripleProduct(const Vector3D& v1, const Vector3D& v2);
	EXPORT float VectorGetDeterminant(const Vector2D& v1, const Vector2D& v2);
	EXPORT void PrintProperties(const Vector2D& vector);
	EXPORT void PrintProperties(const Vector3D& vector);

	// Geometry maths

	template<typename T>
	class EXPORT BaseIterator
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

		void AddIntersect(const Intersect& intersect);
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

		void Translate(const Float2& translation);
		void Translate(const float translateX, const float translateY);
		void Rotate(const float angle);
		void Scale(const float scale);
		void Scale(const float scaleX, const float scaleY);

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

		void Translate(const Float2& translation);
		void Translate(const float translateX, const float translateY);

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
		return round(Pow(10, decimal) * num) / Pow(10, decimal);
	}

	template <typename T>
	void SetArraySize(T** arr, const size_t& curLen, const size_t& newLen)
	{
		T* newArr = new T[newLen];
		memcpy(newArr, *arr, Min(curLen, newLen) * sizeof(**arr));
		delete[] *arr;
		*arr = newArr;
	}
}
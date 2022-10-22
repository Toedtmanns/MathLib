#pragma once
#include "mlPrimitives.hpp"
#include "mlMatrices.hpp"

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
		Line2D(const double x1, const double y1, const double x2, const double y2);
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

	EXPORT double GetLength(const Line2D& line);
	EXPORT double GetDistance(const Float2& p1, const Float2& p2);
	EXPORT double GetDistance(const Float3& p1, const Float3& p2);

	EXPORT double GetSlope(const Line2D& line);
	EXPORT double GetYAxisSection(const Line2D& line);
	EXPORT bool IsParallel(const Line2D& l1, const Line2D& l2);
	EXPORT Float2 GetIntersectPos(const Line2D& l1, const Line2D& l2);
	EXPORT Intersect GetIntersect(const Line2D& l1, const Line2D& l2);
	EXPORT Intersect GetIntersect(const Line2D& l, const Float2& p);

	EXPORT void PrintProperties(const Line2D& l);
	EXPORT void PrintProperties(const Intersect& i);

	double Lerp(const double start, const double end, const double t);
	Float2 Lerp(const Float2& start, const Float2& end, const double t);
	Float2 Lerp(const Line2D& line, const double t);

	// Quaternions

	class EXPORT Quaternion
	{
	public:
		double real, i, j, k;

		Quaternion();
		Quaternion(const double real, const double i, const double j, const double k);
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

	EXPORT Quaternion QuaternionRotation(const double angle, const double iAxis, const double jAxis, const double kAxis);
	EXPORT Quaternion QuaternionRotation(const double angle, const Float3& axis);
	EXPORT void PrintProperties(const Quaternion& quat);

	// Quaternion operations

	EXPORT Mat4 MatrixRotate(const Mat4& mat, const Quaternion& quat);

	// Conversion and checking functions

	EXPORT bool MatrixIsSquare(const MatrixF& mat, const int& dimension = -1);

	// Helper functions

	EXPORT MatrixF Point2Matrix(const Float2& point);
	EXPORT MatrixF Point2Matrix(const Float3& point);
	EXPORT MatrixF TransformF2x2(const Float2& p1, const Float2& p2);

	// Vector math

	class EXPORT Vector2D : public Float2
	{
	public:
		Vector2D();
		Vector2D(const Float2& dir);
		explicit Vector2D(const double x, const double y);
		Vector2D(const Line2D& line);
		Vector2D(const Float2& p1, const Float2& p2);

		const Vector2D& Transform(const Mat2& transformMat);
		const Vector2D& Scale(const double scale);
		const Vector2D& Scale(const double scaleX, const double scaleY);
		const Vector2D& SetScale(const double scale);
		const Vector2D& Rotate(double angle);
		const Vector2D& Rot90R();
		const Vector2D& Rot90L();
		const Vector2D& Normalize();

		Vector2D TripleProduct(const Vector2D& other) const;

		Vector2D operator+(const Vector2D& other) const;
		Vector2D operator-(const Vector2D& other) const;
		double operator*(const Vector2D& other) const;
		Vector2D operator*(const double number) const;
		Vector2D operator*(const Mat2& matrix) const;
		void operator*=(const double number);
		Vector2D operator-() const;

		double GetAngle() const;
		double GetLen() const;
		Mat2x1 GetRowVector() const;
		Mat1x2 GetColVector() const;

		Float2 TransformPoint(Float2 point) const;
	};

	class EXPORT Vector3D : public Float3
	{
	public:
		Vector3D();
		Vector3D(const Float3& dir);
		explicit Vector3D(const double x, const double y, const double z);
		Vector3D(const Vector2D& vec, const double z = 0);
		Vector3D(const Float3& p1, const Float3& p2);

		const Vector3D& Transform(const Mat3& mat);
		const Vector3D& Scale(const double& scale);
		const Vector3D& Scale(const double& scaleX, const double& scaleY, const double& scaleZ);
		const Vector3D& SetScale(const double& len);
		const Vector3D& Rotate(const double& angle, const unsigned int& axis);
		const Vector3D& Normalize();

		Vector3D operator+(const Vector3D& other) const;
		Vector3D operator-(const Vector3D& other) const;
		double operator*(const Vector3D& other) const;
		Vector3D operator*(const double number) const;
		Vector3D operator*(const Mat3& matrix) const;
		void operator*=(const double number);
		Vector3D operator-() const;

		double GetAngle(const int axis) const;
		double GetLen() const;
		Mat3x1 GetRowVector() const;
		Mat1x3 GetColVector() const;

		Float3 TransformPoint(Float3 point);
	};

	EXPORT Vector2D operator*(const Mat2& matrix, const Vector2D& vector);
	EXPORT Vector3D operator*(const Mat3& matrix, const Vector3D& vector);

	EXPORT Vector3D GetRelativeVec(const Vector3D& vec1, const Vector3D& vec2);
	EXPORT double VectorGetAngleDifference(const Vector2D& v1, const Vector2D& v2);
	EXPORT double VectorDotProduct(const Vector2D& v1, const Vector2D& v2);
	EXPORT Vector3D VectorCrossProduct(const Vector2D& v1, const Vector2D& v2);
	EXPORT Vector3D VectorCrossProduct(const Vector3D& v1, const Vector3D& v2);
	EXPORT Vector2D VectorTripleProduct(const Vector2D& v1, const Vector2D& v2);
	EXPORT Vector3D VectorTripleProduct(const Vector3D& v1, const Vector3D& v2);
	EXPORT double VectorGetDeterminant(const Vector2D& v1, const Vector2D& v2);
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
		void Translate(const double translateX, const double translateY);
		void Rotate(const double angle);
		void Scale(const double scale);
		void Scale(const double scaleX, const double scaleY);

		Float2 GetCenter() const;
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
		double radius;

		Circle2D();
		Circle2D(const Float2& position, const double radius = 0.5);

		void Translate(const Float2& translation);
		void Translate(const double translateX, const double translateY);

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
		Rectangle2D(const Float2& position, const double rotation, const Float2& scale);

		void operator=(Rectangle2D&& other) noexcept;
		void operator=(const Rectangle2D& other);
	};

	EXPORT bool Contains(const Rectangle2D& rect, const double rotation, const Float2& point);
	EXPORT double* ProjectTo1D(const Polygon2D& polygon, const Vector2D& viewDir);
	EXPORT double* ProjectTo1D(const Float2* pointArr, const size_t pointCount, const Vector2D& viewDir);

	// Utility
	EXPORT double MinFromArray(const double* arr, const size_t length);
	EXPORT double MaxFromArray(const double* arr, const size_t length);
	EXPORT Line2D Vector2Line(Vector2D vector, Float2 pos);
	EXPORT int RandInt(int min, int max);
	EXPORT inline double RoundFTo(const double num, const int decimal)
	{
		return round(pow(10, decimal) * num) / pow(10, decimal);
	}
	EXPORT MatrixF Vec2Mat(const Vector3D& vec, const int front = 0);
	EXPORT Vector3D Mat2Vec(const MatrixF& mat, const int front = 0);

	EXPORT bool SetArraySize(void** array, const size_t currSize, const size_t newSize);
}
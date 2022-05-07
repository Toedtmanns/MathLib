#pragma once
#define MATHLIB

#include <cmath>

#define PI 3.14159265

#ifndef MATHLIB_STATIC
#ifndef EXPORT
#ifdef MATHLIB_EXPORTS
#define EXPORT __declspec(dllexport)
#else
#define EXPORT __declspec(dllimport)
#endif // MATHLIB_EXPORTS
#endif // EXPORT

#else
#ifndef EXPORT
#define EXPORT
#endif // EXPORT
#endif // MATHLIB_STATIC

namespace MathLib
{
	// Generally useful functions

	template<typename T>
	T Max(T t1, T t2)
	{
		if (t1 > t2)
			return t1;
		return t2;
	}

	template<typename T>
	T Min(T t1, T t2)
	{
		if (t1 < t2)
			return t1;
		return t2;
	}

	// Primitive types

	class Int4;
	class Int3;
	class Int2;

	class EXPORT Float4
	{
	public:
		union
		{
			struct
			{
				double x;
				double y;
				double z;
				double w;
			};
			struct
			{
				double r;
				double g;
				double b;
				double a;
			};
			double content[4];
		};

		Float4();
		Float4(const double& x, const double& y, const double& z, const double& w);
		Float4(const Int4& other);

		void operator=(const Int4& other);
		double& operator[](const unsigned int& index);
		const double& operator[](const unsigned int& index) const;
	};

	class EXPORT Float3
	{
	public:
		union
		{
			struct
			{
				double x;
				double y;
				double z;
			};
			struct
			{
				double r;
				double g;
				double b;
			};
			double content[3];
		};

		Float3();
		Float3(const double& x, const double& y, const double& z);
		Float3(const Int3& other);

		void operator=(const Int3& other);
		double& operator[](const unsigned int& index);
		const double& operator[](const unsigned int& index) const;
	};

	class EXPORT Float2
	{
	public:
		double x;
		double y;

		Float2();
		Float2(const double& x, const double& y);
		Float2(const Int2& other);

		void operator=(const Int2& other);
		double& operator[](const unsigned int& index);
		const double& operator[](const unsigned int& index) const;
	};

	class EXPORT Int4
	{
	public:
		union
		{
			struct
			{
				int x;
				int y;
				int z;
				int w;
			};
			struct
			{
				int r;
				int g;
				int b;
				int a;
			};
			int content[4];
		};

		Int4();
		Int4(const int& x, const int& y, const int& z, const int& w);
		Int4(const Float4& other);

		void operator=(const Float4& other);
		int& operator[](const unsigned int& index);
		const int& operator[](const unsigned int& index) const;
	};

	class EXPORT Int3
	{
	public:
		union
		{
			struct
			{
				int x;
				int y;
				int z;
			};
			struct
			{
				int r;
				int g;
				int b;
			};
			int content[3];
		};

		Int3();
		Int3(const int& x, const int& y, const int& z);
		Int3(const Float3& other);

		void operator=(const Float3& other);
		int& operator[](const unsigned int& index);
		const int& operator[](const unsigned int& index) const;
	};

	class EXPORT Int2
	{
	public:
		int x;
		int y;

		Int2();
		Int2(const int& x, const int& y);
		Int2(const Float2& other);

		void operator=(const Float2& other);
		int& operator[](const unsigned int& index);
		const int& operator[](const unsigned int& index) const;
	};

	EXPORT bool operator==(const Float2& f1, const Float2& f2);
	EXPORT bool operator!=(const Float2& f1, const Float2& f2);
	EXPORT bool operator==(const Float3& f1, const Float3& f2);
	EXPORT bool operator!=(const Float3& f1, const Float3& f2);
	EXPORT bool operator==(const Float4& f1, const Float4& f2);
	EXPORT bool operator!=(const Float4& f1, const Float4& f2);
	EXPORT Float2 operator+(const Float2& f1, const Float2& f2);
	EXPORT Float3 operator+(const Float3& f1, const Float3& f2);
	EXPORT Float4 operator+(const Float4& f1, const Float4& f2);
	EXPORT Float2 operator-(const Float2& f1, const Float2& f2);
	EXPORT Float3 operator-(const Float3& f1, const Float3& f2);
	EXPORT Float4 operator-(const Float4& f1, const Float4& f2);
	EXPORT void operator+=(Float2& f1, const Float2& f2);
	EXPORT void operator+=(Float3& f1, const Float3& f2);
	EXPORT void operator+=(Float4& f1, const Float4& f2);
	EXPORT void operator-=(Float2& f1, const Float2& f2);
	EXPORT void operator-=(Float3& f1, const Float3& f2);
	EXPORT void operator-=(Float4& f1, const Float4& f2);
	EXPORT Float2 operator*(const Float2& f1, const Float2& f2);
	EXPORT Float3 operator*(const Float3& f1, const Float3& f2);
	EXPORT Float4 operator*(const Float4& f1, const Float4& f2);
	EXPORT Float2 operator*(const Float2& f1, const double& f2);
	EXPORT Float3 operator*(const Float3& f1, const double& f2);
	EXPORT Float4 operator*(const Float4& f1, const double& f2);

	EXPORT bool operator==(const Int2& i1, const Int2& i2);
	EXPORT bool operator!=(const Int2& i1, const Int2& i2);
	EXPORT bool operator==(const Int3& i1, const Int3& i2);
	EXPORT bool operator!=(const Int3& i1, const Int3& i2);
	EXPORT bool operator==(const Int4& i1, const Int4& i2);
	EXPORT bool operator!=(const Int4& i1, const Int4& i2);
	EXPORT Int2 operator+(const Int2& i1, const Int2& i2);
	EXPORT Int3 operator+(const Int3& i1, const Int3& i2);
	EXPORT Int4 operator+(const Int4& i1, const Int4& i2);
	EXPORT Int2 operator-(const Int2& i1, const Int2& i2);
	EXPORT Int3 operator-(const Int3& i1, const Int3& i2);
	EXPORT Int4 operator-(const Int4& i1, const Int4& i2);
	EXPORT void operator+=(Int2& i1, const Int2& i2);
	EXPORT void operator+=(Int3& i1, const Int3& i2);
	EXPORT void operator+=(Int4& i1, const Int4& i2);
	EXPORT void operator-=(Int2& i1, const Int2& i2);
	EXPORT void operator-=(Int3& i1, const Int3& i2);
	EXPORT void operator-=(Int4& i1, const Int4& i2);
	EXPORT Int2 operator*(const Int2& f1, const Int2& f2);
	EXPORT Int3 operator*(const Int3& f1, const Int3& f2);
	EXPORT Int4 operator*(const Int4& f1, const Int4& f2);
	EXPORT Int2 operator*(const Int2& f1, const int& f2);
	EXPORT Int3 operator*(const Int3& f1, const int& f2);
	EXPORT Int4 operator*(const Int4& f1, const int& f2);

	EXPORT void PrintProperties(const Float2& p);
	EXPORT void PrintProperties(const Float3& p);

	// Lines

	struct EXPORT Line2D
	{
		Float2 p1;
		Float2 p2;
		double normal;

		Line2D();
		Line2D(const double& x1, const double& y1, const double& x2, const double& y2);
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
		Intersect(const Float2& pos, const bool& intersecting, const bool& collinear);
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

	double Lerp(const double& start, const double& end, const double& t);
	Float2 Lerp(const Float2& start, const Float2& end, const double& t);
	Float2 Lerp(const Line2D& line, const double& t);

	// Complex maths

	class EXPORT imaginaryBase
	{
	public:
		double num;

		imaginaryBase();
		imaginaryBase(const double& num);

		template <typename T>
		imaginaryBase operator*(const T& other)
		{
			return imaginaryBase(num * other);
		}
		imaginaryBase operator+(const imaginaryBase& other) const;
		imaginaryBase operator-(const imaginaryBase& other) const;
		double operator*(const imaginaryBase& other) const;
		double operator/(const imaginaryBase& other) const;
		bool operator==(const imaginaryBase& other) const;
		bool operator!=(const imaginaryBase& other) const;
		bool operator>=(const imaginaryBase& other) const;
		bool operator<=(const imaginaryBase& other) const;
		bool operator>(const imaginaryBase& other) const;
		bool operator<(const imaginaryBase& other) const;
		imaginaryBase& operator++();
		imaginaryBase& operator--();
		imaginaryBase& operator++(int);
		imaginaryBase& operator--(int);
	};

	template <typename T>
	imaginaryBase operator*(const T& num1, const imaginaryBase& num2)
	{
		return imaginaryBase(num1 * num2.num);
	}
	EXPORT imaginaryBase operator-(const imaginaryBase& num);

	class imagJ;
	class imagK;

	class EXPORT imagI : public imaginaryBase
	{
	public:
		imagI();
		imagI(const double& num);
		imagI(const imaginaryBase& base);

		double operator*(const imagI& other) const;
		imagK operator*(const imagJ& jNum) const;
		imagJ operator*(const imagK& jNum) const;
	};
	class EXPORT imagJ : public imaginaryBase
	{
	public:
		imagJ();
		imagJ(const double& num);
		imagJ(const imaginaryBase& base);

		double operator*(const imagJ& other) const;
		imagI operator*(const imagK& jNum) const;
		imagK operator*(const imagI& jNum) const;
	};
	class EXPORT imagK : public imaginaryBase
	{
	public:
		imagK();
		imagK(const double& num);
		imagK(const imaginaryBase& base);

		double operator*(const imagK& other) const;
		imagJ operator*(const imagI& jNum) const;
		imagI operator*(const imagJ& jNum) const;
	};

	// Quaternions

	class EXPORT Quaternion
	{
	public:
		double real;
		imagI i;
		imagJ j;
		imagK k;

		Quaternion();
		Quaternion(const double& real, const imagI& i, const imagJ& j, const imagK& k);
		Quaternion(const Float3& point);

		Quaternion operator+(const Quaternion& other) const;
		Quaternion operator*(const Quaternion& other) const;
		Quaternion RotateQuaternion(const Quaternion& quat) const;
		Float3 RotatePoint(const Float3& point) const;

		Quaternion GetInverse() const;
		Float3 GetPoint() const;
	};

	EXPORT Quaternion QuaternionRotation(const double& angle, const double& iAxis, const double& jAxis, const double& kAxis);
	EXPORT Quaternion QuaternionRotation(const double& angle, const Float3& axis);
	EXPORT void PrintProperties(const Quaternion& quat);

	// Matrix maths

	class EXPORT MatrixI
	{
		unsigned int m_Rows, m_Columns;
		int** m_Matrix;

	public:
		MatrixI(MatrixI&& other) noexcept;
		MatrixI(const MatrixI& other);
		MatrixI(const unsigned int& dim);
		MatrixI(const unsigned int& columns, const unsigned int& rows);
		MatrixI(const unsigned int& columns, const unsigned int& rows, const int** const columnArray);

		void SetNum(const unsigned int& column, const unsigned int& row, const int& value);
		void SetColumn(const unsigned int& column, const int* const content);
		void SetMatrix(const int** const matrix);

		const int& GetNum(const unsigned int& column, const unsigned int& row) const;
		int* GetColumn(const unsigned int& column) const;
		int* GetArray() const;
		int** GetMatrix() const;
		const unsigned int& GetRowCount() const;
		const unsigned int& GetColumnCount() const;

		void operator=(MatrixI&& other) noexcept;
		void operator=(const MatrixI& other);
		MatrixI operator+(const int& value) const;
		MatrixI operator-(const int& value) const;
		MatrixI operator*(const int& value) const;
		MatrixI operator/(const int& value) const;
		MatrixI operator+(const MatrixI& other) const;
		MatrixI operator-(const MatrixI& other) const;
		MatrixI operator*(const MatrixI& other) const;
		int* operator[](const unsigned int& column);
		const int* operator[](const unsigned int& column) const;

		~MatrixI();
	};
	class EXPORT MatrixF
	{
		double** m_Matrix;
		unsigned int m_Rows, m_Columns;

	public:
		MatrixF(MatrixF&& other) noexcept;
		MatrixF(const MatrixF& other);
		MatrixF(const unsigned int& dim);
		MatrixF(const unsigned int& columns, const unsigned int& rows);
		MatrixF(const unsigned int& columns, const unsigned int& rows, const double** const columnArray);

		void SetNum(const unsigned int& column, const unsigned int& row, const double& value);
		void SetColumn(const unsigned int& column, const double* const content);
		void SetMatrix(const double** const matrix);

		const double& GetNum(const unsigned int& column, const unsigned int& row) const;
		double* GetColumn(const unsigned int& column) const;
		double* GetArray() const;
		double** GetMatrix() const;
		const unsigned int& GetRowCount() const;
		const unsigned int& GetColumnCount() const;

		void operator=(MatrixF&& other) noexcept;
		void operator=(const MatrixF& other);
		MatrixF operator+(const double& value) const;
		MatrixF operator-(const double& value) const;
		MatrixF operator*(const double& value) const;
		MatrixF operator/(const double& value) const;
		MatrixF operator+(const MatrixF& other) const;
		MatrixF operator-(const MatrixF& other) const;
		MatrixF operator*(const MatrixF& other) const;
		double* operator[](const unsigned int& column);
		const double* operator[](const unsigned int& column) const;

		~MatrixF();
	};

	EXPORT void PrintContent(const MatrixI& mat);
	EXPORT void PrintProperties(const MatrixI& mat);
	EXPORT int MatrixGetDet(const MatrixI& mat);
	EXPORT MatrixI MatrixOfMinors(const MatrixI& mat);
	EXPORT MatrixI MatrixOfCofactors(const MatrixI& mat);
	EXPORT MatrixI MatrixTranspose(const MatrixI& mat);
	EXPORT MatrixI MatrixAdjugate(const MatrixI& mat);
	EXPORT MatrixF MatrixInverse(const MatrixI& mat);

	EXPORT void PrintContent(const MatrixF& mat);
	EXPORT void PrintProperties(const MatrixF& mat);
	EXPORT double MatrixGetDet(const MatrixF& mat);
	EXPORT MatrixF MatrixOfMinors(const MatrixF& mat);
	EXPORT MatrixF MatrixOfCofactors(const MatrixF& mat);
	EXPORT MatrixF MatrixTranspose(const MatrixF& mat);
	EXPORT MatrixF MatrixAdjugate(const MatrixF& mat);
	EXPORT MatrixF MatrixInverse(const MatrixF& mat);

	// Quaternion operations

	EXPORT MatrixF MatrixRotate(const MatrixF& mat, const Quaternion& quat);

	// Conversion and checking functions

	EXPORT MatrixI MatrixF2I(const MatrixF& mat);
	EXPORT MatrixF MatrixI2F(const MatrixI& mat);
	EXPORT bool MatrixIsSquare(const MatrixI& mat, const int& dimension = -1);
	EXPORT bool MatrixIsSquare(const MatrixF& mat, const int& dimension = -1);

	// Helper functions

	EXPORT MatrixI Point2Matrix(const Int2& point);
	EXPORT MatrixF Point2Matrix(const Float2& point);
	EXPORT MatrixI Point2Matrix(const Int3& point);
	EXPORT MatrixF Point2Matrix(const Float3& point);
	EXPORT MatrixF TransformF2x2(const Float2& p1, const Float2& p2);
	EXPORT MatrixI TransformI2x2(const Int2& p1, const Int2& p2);

	// Vector math

	class EXPORT Vector2D : public Float2
	{
	public:
		Vector2D();
		Vector2D(const Float2& dir);
		Vector2D(const double& x, const double& y);
		Vector2D(const Line2D& line);
		Vector2D(const Float2& p1, const Float2& p2);

		const Vector2D& Transform(const MatrixF& transformMat);
		const Vector2D& Scale(const double& scale);
		const Vector2D& Scale(const double& scaleX, const double& scaleY);
		const Vector2D& SetScale(const double& scale);
		const Vector2D& Rotate(double angle);
		const Vector2D& Rot90R();
		const Vector2D& Rot90L();

		Vector2D TripleProduct(const Vector2D& other) const;

		Vector2D operator+(const Vector2D& other) const;
		Vector2D operator-(const Vector2D& other) const;
		double operator*(const Vector2D& other) const;
		Vector2D operator*(const double& number) const;
		Vector2D operator*(const MatrixF& matrix) const;
		void operator*=(const double& number);
		Vector2D operator-() const;

		double GetAngle() const;
		double GetLen() const;
		MatrixF GetRowVector() const;
		MatrixF GetColVector() const;

		Float2 TransformPoint(Float2 point) const;
	};

	class EXPORT Vector3D : public Float3
	{
	public:
		Vector3D();
		Vector3D(const Float3& dir);
		Vector3D(const double& x, const double& y, const double& z);
		Vector3D(const Vector2D& vec, const double& z = 0);
		Vector3D(const Float3& p1, const Float3& p2);

		const Vector3D& Transform(const MatrixF& mat);
		const Vector3D& Scale(const double& scale);
		const Vector3D& Scale(const double& scaleX, const double& scaleY, const double& scaleZ);
		const Vector3D& SetScale(const double& len);
		const Vector3D& Rotate(const double& angle, const unsigned int& axis);

		Vector3D operator+(const Vector3D& other) const;
		Vector3D operator-(const Vector3D& other) const;
		double operator*(const Vector3D& other) const;
		Vector3D operator*(const double& number) const;
		Vector3D operator*(const MatrixF& matrix) const;
		void operator*=(const double& number);
		Vector3D operator-() const;

		double GetAngle(const int& axis) const;
		double GetLen() const;
		MatrixF GetRowVector() const;
		MatrixF GetColVector() const;

		Float3 TransformPoint(Float3 point);
	};

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

	EXPORT Float2 operator*(const Float2& point, const MatrixF& matrix);
	EXPORT Float2 Transform(Float2 point, Float2 origin, MatrixF* transform);
	EXPORT Line2D Transform(Line2D line, Float2 origin, MatrixF* transform);
	EXPORT Float3 RotatePoint(Float3 point, const Quaternion& quat);

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
		ValType& operator[](const unsigned int& index)
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
		unsigned int m_IntersectArrLength;
		unsigned int m_IntersectCount;
		Intersect* m_IntersectArray;
	public:
		using Iterator = BaseIterator<GeoCollision>;
		using ValueType = Intersect;
		bool isColliding;

		GeoCollision();
		GeoCollision(GeoCollision&& other) noexcept;
		GeoCollision(const GeoCollision& other);
		GeoCollision(const Intersect* intersectArr, const unsigned int& intersectArrLength);

		void AddIntersect(const Intersect& intersect);
		unsigned int GetIntersectCount() const;

		void operator=(GeoCollision&& other) noexcept;
		void operator=(const GeoCollision& other);
		Intersect& operator[](const unsigned int& index);

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
		unsigned int m_NumCorners;
		Float2* m_CornerArr;

		Polygon2D();
		Polygon2D(Polygon2D&& other) noexcept;
		Polygon2D(const Polygon2D& other);
		Polygon2D(const unsigned int& corners);
		Polygon2D(const Float2* const pointArr, const unsigned int& corners);

		void Translate(const Float2& translation);
		void Translate(const double& translateX, const double& translateY);
		void Rotate(const double& angle);
		void Scale(const double& scale);
		void Scale(const double& scaleX, const double& scaleY);

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
		Circle2D(const Float2& position, const double& radius = 0.5);

		void Translate(const Float2& translation);
		void Translate(const double& translateX, const double& translateY);

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

		bool CollidesWith(const Triangle2D& other) const;
		GeoCollision GetCollision(const Triangle2D& other) const;
	};

	class EXPORT Rectangle2D : public Polygon2D
	{
	public:
		Rectangle2D();
		Rectangle2D(Rectangle2D&& other) noexcept;
		Rectangle2D(const Rectangle2D& other);
		Rectangle2D(const Float2& p1, const Float2& p2, const Float2& p3, const Float2& p4);
		Rectangle2D(const Float2* const pointArr);
		Rectangle2D(const Float2 position, const double& rotation, const Float2 scale);

		void operator=(Rectangle2D&& other) noexcept;
		void operator=(const Rectangle2D& other);

		bool CollidesWith(const Rectangle2D& other) const;
		GeoCollision GetCollision(const Rectangle2D& other) const;
	};

	EXPORT bool Contains(const Rectangle2D& rect, const double& rotation, const Float2& point);
	EXPORT double* ProjectTo1D(const Polygon2D& polygon, const Vector2D& viewDir);
	EXPORT double* ProjectTo1D(const Float2* pointArr, const unsigned int& pointCount, const Vector2D& viewDir);

	// Utility

	EXPORT double Deg2Rad(double deg);
	EXPORT double Rad2Deg(double rad);
	EXPORT double MinFromArray(const double* const arr, const unsigned int& length);
	EXPORT double MaxFromArray(const double* const arr, const unsigned int& length);
	EXPORT Line2D Vector2Line(Vector2D vector, Float2 pos);
	EXPORT int RandInt(int min, int max);
	EXPORT inline double RoundFTo(const double& num, const int& decimal)
	{
		return round(pow(10, decimal) * num) / pow(10, decimal);
	}
	EXPORT MatrixF Vec2Mat(const Vector3D& vec, const int& front = 0);
	EXPORT Vector3D Mat2Vec(const MatrixF& mat, const int& front = 0);

	EXPORT bool SetArraySize(void** array, const unsigned int& currSize, const unsigned int& newSize);
}
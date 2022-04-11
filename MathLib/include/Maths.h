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

	namespace Primitives
	{
		union Int4;
		union Int3;
		union Int2;

		union EXPORT Float4
		{
		public:
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

			Float4();
			Float4(const double& x, const double& y, const double& z, const double& w);
			Float4(const Int4& other);

			void operator=(const Int4& other);
			double& operator[](const unsigned int& index);
			const double& operator[](const unsigned int& index) const;
		};

		union EXPORT Float3
		{
		public:
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

			Float3();
			Float3(const double& x, const double& y, const double& z);
			Float3(const Int3& other);

			void operator=(const Int3& other);
			double& operator[](const unsigned int& index);
			const double& operator[](const unsigned int& index) const;
		};

		union EXPORT Float2
		{
		public:
			struct
			{
				double x;
				double y;
			};
			double content[2];

			Float2();
			Float2(const double& x, const double& y);
			Float2(const Int2& other);

			void operator=(const Int2& other);
			double& operator[](const unsigned int& index);
			const double& operator[](const unsigned int& index) const;
		};

		union EXPORT Int4
		{
		public:
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

			Int4();
			Int4(const int& x, const int& y, const int& z, const int& w);
			Int4(const Float4& other);

			void operator=(const Float4& other);
			int& operator[](const unsigned int& index);
			const int& operator[](const unsigned int& index) const;
		};

		union EXPORT Int3
		{
		public:
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

			Int3();
			Int3(const int& x, const int& y, const int& z);
			Int3(const Float3& other);

			void operator=(const Float3& other);
			int& operator[](const unsigned int& index);
			const int& operator[](const unsigned int& index) const;
		};

		union EXPORT Int2
		{
		public:
			struct
			{
				int x;
				int y;
			};
			int content[2];

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
	}

	double Lerp(const double& start, const double& end, const double& t);
	Primitives::Float2 Lerp(const Primitives::Float2& start, const Primitives::Float2& end, const double& t);
	Primitives::Float2 Lerp(const Primitives::Line2D& line, const double& t);

	namespace Complex
	{
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
			Quaternion(const Primitives::Float3& point);

			Quaternion operator+(const Quaternion& other) const;
			Quaternion operator*(const Quaternion& other) const;
			Quaternion RotateQuaternion(const Quaternion& quat) const;
			Primitives::Float3 RotatePoint(const Primitives::Float3& point) const;

			Quaternion GetInverse() const;
			Primitives::Float3 GetPoint() const;
		};

		EXPORT Quaternion QuaternionRotation(const double& angle, const double& iAxis, const double& jAxis, const double& kAxis);
		EXPORT Quaternion QuaternionRotation(const double& angle, const Primitives::Float3& axis);
		EXPORT void PrintProperties(const Quaternion& quat);
	}

	namespace Matrices
	{
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

		EXPORT MatrixF MatrixRotate(const MatrixF& mat, const Complex::Quaternion& quat);

		// Conversion and checking functions

		EXPORT MatrixI MatrixF2I(const MatrixF& mat);
		EXPORT MatrixF MatrixI2F(const MatrixI& mat);
		EXPORT bool MatrixIsSquare(const MatrixI& mat, const int& dimension = -1);
		EXPORT bool MatrixIsSquare(const MatrixF& mat, const int& dimension = -1);

		// Helper functions

		EXPORT MatrixI Point2Matrix(const Primitives::Int2& point);
		EXPORT MatrixF Point2Matrix(const Primitives::Float2& point);
		EXPORT MatrixI Point2Matrix(const Primitives::Int3& point);
		EXPORT MatrixF Point2Matrix(const Primitives::Float3& point);
		EXPORT MatrixI TransformI2x2(const Primitives::Int2& p1, const Primitives::Int2& p2);
		EXPORT MatrixF TransformF2x2(const Primitives::Float2& p1, const Primitives::Float2& p2);
	}

	namespace Vectors
	{
		class EXPORT Vector2D
		{
		public:
			Primitives::Float2 direction;

			Vector2D();
			Vector2D(const Primitives::Float2& dir);
			Vector2D(const double& angle, const double& length);
			Vector2D(const Primitives::Line2D& line);

			void Transform(const Matrices::MatrixF& transformMat);
			void Scale(const double& scale);
			void Scale(const double& scaleX, const double& scaleY);
			void SetScale(const double& scale);
			void Rotate(double angle);

			Vector2D operator+(const Vector2D& other) const;
			Vector2D operator-(const Vector2D& other) const;
			double operator*(const Vector2D& other) const;

			double GetAngle() const;
			double GetLen() const;
			Matrices::MatrixF GetRowVector() const;
			Matrices::MatrixF GetColVector() const;

			Primitives::Float2 TransformPoint(Primitives::Float2 point) const;
		};

		class EXPORT Vector3D
		{
		public:
			Primitives::Float3 direction;

			Vector3D();
			Vector3D(const Primitives::Float3& dir);

			void Transform(const Matrices::MatrixF& mat);
			void Scale(const double& scale);
			void Scale(const double& scaleX, const double& scaleY, const double& scaleZ);
			void Rotate(const double& angle, const unsigned int& axis);
			void SetLen(const double& len);

			Vector3D operator+(const Vector3D& other) const;
			Vector3D operator-(const Vector3D& other) const;
			double operator*(const Vector3D& other) const;

			double GetAngle(const int& axis) const;
			//double GetAngle(const Vector3D& normal) const;
			double GetLen() const;
			Matrices::MatrixF GetRowVector() const;
			Matrices::MatrixF GetColVector() const;

			Primitives::Float3 TransformPoint(Primitives::Float3 point);
		};

		EXPORT Vector3D GetRelativeVec(const Vector3D& vec1, const Vector3D& vec2);
		EXPORT double VectorGetAngleDifference(const Vector2D& v1, const Vector2D& v2);
		EXPORT double VectorDotProduct(const Vector2D& v1, const Vector2D& v2);
		EXPORT Vector3D VectorCrossProduct(const Vector2D& v1, const Vector2D& v2);
		EXPORT double VectorGetDeterminant(const Vector2D& v1, const Vector2D& v2);
		EXPORT void PrintProperties(const Vector2D& vector);
		EXPORT void PrintProperties(const Vector3D& vector);
	}

	namespace Primitives
	{
		EXPORT Float2 operator*(const Float2& point, const Matrices::MatrixF& matrix);
		EXPORT Float2 Transform(Float2 point, Float2 origin, Matrices::MatrixF* transform);
		EXPORT Line2D Transform(Line2D line, Float2 origin, Matrices::MatrixF* transform);
		EXPORT Float3 RotatePoint(Float3 point, const Complex::Quaternion& quat);
	}

	namespace Geometry
	{
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
			Primitives::Intersect* m_IntersectArray;
		public:
			using Iterator = BaseIterator<GeoCollision>;
			using ValueType = Primitives::Intersect;
			bool isColliding;

			GeoCollision();
			GeoCollision(GeoCollision&& other) noexcept;
			GeoCollision(const GeoCollision& other);
			GeoCollision(const Primitives::Intersect* intersectArr, const unsigned int& intersectArrLength);

			void AddIntersect(const Primitives::Intersect& intersect);
			unsigned int GetIntersectCount() const;

			void operator=(GeoCollision&& other) noexcept;
			void operator=(const GeoCollision& other);
			Primitives::Intersect& operator[](const unsigned int& index);

			Iterator begin();
			Iterator end();

			~GeoCollision();
		};

		class EXPORT Polygon2D
		{
		public:
			using Iterator = BaseIterator<Polygon2D>;
			using ValueType = Primitives::Float2;
			unsigned int m_Corners;
			Primitives::Float2* m_PointArr;

			Polygon2D();
			Polygon2D(Polygon2D&& other) noexcept;
			Polygon2D(const Polygon2D& other);
			Polygon2D(const unsigned int& corners);
			Polygon2D(const Primitives::Float2* const pointArr, const unsigned int& corners);

			void Translate(const Primitives::Float2& translation);
			void Translate(const double& translateX, const double& translateY);
			void Rotate(const double& angle);
			void Scale(const double& scale);
			void Scale(const double& scaleX, const double& scaleY);

			Primitives::Float2 GetCenter() const;

			void operator=(Polygon2D&& other) noexcept;
			void operator=(const Polygon2D& other);

			Iterator begin();
			Iterator end();

			~Polygon2D();
		};

		class EXPORT Triangle2D : public Polygon2D
		{
		public:
			Triangle2D();
			Triangle2D(Triangle2D&& other) noexcept;
			Triangle2D(const Triangle2D& other);
			Triangle2D(const Primitives::Float2& p1, const Primitives::Float2& p2, const Primitives::Float2& p3);
			Triangle2D(const Primitives::Float2* const pointArr);

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
			Rectangle2D(const Primitives::Float2& p1, const Primitives::Float2& p2, const Primitives::Float2& p3, const Primitives::Float2& p4);
			Rectangle2D(const Primitives::Float2* const pointArr);
			Rectangle2D(const Primitives::Float2 position, const double& rotation, const Primitives::Float2 scale);

			void operator=(Rectangle2D&& other) noexcept;
			void operator=(const Rectangle2D& other);

			bool CollidesWith(const Rectangle2D& other) const;
			GeoCollision GetCollision(const Rectangle2D& other) const;
		};

		EXPORT bool Contains(const Rectangle2D& rect, const double& rotation, const Primitives::Float2& point);
		EXPORT double* ProjectTo1D(const Polygon2D& polygon, const Vectors::Vector2D& viewDir);
		EXPORT double* ProjectTo1D(const Primitives::Float2* pointArr, const unsigned int& pointCount, const Vectors::Vector2D& viewDir);
	}

	namespace Utility
	{
		EXPORT double Deg2Rad(double deg);
		EXPORT double Rad2Deg(double rad);
		EXPORT double MinFromArray(const double* const arr, const unsigned int& length);
		EXPORT double MaxFromArray(const double* const arr, const unsigned int& length);
		EXPORT Primitives::Line2D Vector2Line(Vectors::Vector2D vector, Primitives::Float2 pos);
		EXPORT Vectors::Vector2D Line2Vector(Primitives::Line2D line);
		EXPORT Vectors::Vector2D Line2Vector(Primitives::Float2 p1, Primitives::Float2 p2);
		EXPORT int RandInt(int min, int max);
		EXPORT inline double RoundFTo(const double& num, const int& decimal)
		{
			return round(pow(10, decimal) * num) / pow(10, decimal);
		}
		EXPORT Matrices::MatrixF Vec2Mat(const Vectors::Vector3D& vec, const int& front = 0);
		EXPORT Vectors::Vector3D Mat2Vec(const Matrices::MatrixF& mat, const int& front = 0);

		EXPORT bool SetArraySize(void** array, const unsigned int& currSize, const unsigned int& newSize);
	}
}
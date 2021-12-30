#pragma once
#define MATHLIB

#include <cmath>

#ifndef ARRAY_SIZE 
#define ARRAY_SIZE(x) (sizeof((x)) / sizeof((x)[0])
#endif

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
	namespace Primitives
	{
		union EXPORT Float4
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

			Float4()
			{
				x = 0;
				y = 0;
				z = 0;
				w = 0;
			}
			Float4(const double& x, const double& y, const double& z, const double& w)
			{
				this->x = x;
				this->y = y;
				this->z = z;
				this->w = w;
			}
		};

		union EXPORT Float3
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

			Float3()
			{
				x = 0;
				y = 0;
				z = 0;
			}
			Float3(const double& x, const double& y, const double& z)
			{
				this->x = x;
				this->y = y;
				this->z = z;
			}
		};

		struct EXPORT Float2
		{
			double x;
			double y;

			Float2()
			{
				x = 0;
				y = 0;
			}
			Float2(double x, double y)
			{
				this->x = x;
				this->y = y;
			}
		};

		union EXPORT Int4
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

			Int4()
			{
				x = 0;
				y = 0;
				z = 0;
				w = 0;
			}
			Int4(const int& x, const int& y, const int& z, const int& w)
			{
				this->x = x;
				this->y = y;
				this->z = z;
				this->w = w;
			}
		};

		union EXPORT Int3
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

			Int3()
			{
				x = 0;
				y = 0;
				z = 0;
			}
			Int3(const int& x, const int& y, const int& z)
			{
				this->x = x;
				this->y = y;
				this->z = z;
			}
		};

		struct EXPORT Int2
		{
			int x;
			int y;

			Int2()
			{
				x = 0;
				y = 0;
			}
			Int2(const int& x, const int& y)
			{
				this->x = x;
				this->y = y;
			}
		};

		struct EXPORT Line2D
		{
			Float2 p1;
			Float2 p2;
			float normal;

			Line2D();
			Line2D(double x1, double y1, double x2, double y2);
			Line2D(Float2 p1, Float2 p2);
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

		EXPORT void PrintProperties(Float2 p);
		EXPORT void PrintProperties(Float3 p);
		EXPORT void PrintProperties(Line2D l);
	}

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
			int** m_Matrix;
			unsigned int m_Rows, m_Columns;

		public:
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

			void operator=(const MatrixI& other);
			MatrixI operator+(const int& value) const;
			MatrixI operator-(const int& value) const;
			MatrixI operator*(const int& value) const;
			MatrixI operator/(const int& value) const;
			MatrixI operator+(const MatrixI& other) const;
			MatrixI operator-(const MatrixI& other) const;
			MatrixI operator*(const MatrixI& other) const;

			~MatrixI();
		};
		class EXPORT MatrixF
		{
			double** m_Matrix;
			int m_Rows, m_Columns;

		public:
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

			void operator=(const MatrixF& other);
			MatrixF operator+(const double& value) const;
			MatrixF operator-(const double& value) const;
			MatrixF operator*(const double& value) const;
			MatrixF operator/(const double& value) const;
			MatrixF operator+(const MatrixF& other) const;
			MatrixF operator-(const MatrixF& other) const;
			MatrixF operator*(const MatrixF& other) const;

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

		EXPORT MatrixI TransformI2x2(const Primitives::Int2& p1, const Primitives::Int2& p2);
		EXPORT MatrixF TransformF2x2(const Primitives::Float2& p1, const Primitives::Float2& p2);
	}

	namespace Vectors
	{
		class EXPORT Vector3D
		{
		public:
			Primitives::Float3 direction;

			Vector3D();
			Vector3D(Primitives::Float3 dir);

			Vector3D Transform(const Matrices::MatrixF& mat);
			Vector3D Scale(double s);
			Vector3D Scale(double sX, double sY, double sZ);
			Vector3D Rotate(double angle, int axis);
			Vector3D Rotate(double angle, Vector3D* normal);
			Vector3D SetLen(double len);

			double GetAngle(int axis) const;
			double GetAngle(Vector3D* normal) const;
			double GetLen() const;
		};

		/// <summary>
		/// A class which handles vectors and vector calculations
		/// </summary>
		class EXPORT Vector2D
		{
		public:
			Primitives::Float2 direction;

			Vector2D();
			Vector2D(Primitives::Float2 dir);
			Vector2D(double angle, double length);

			Vector2D Transform(const Matrices::MatrixF& mat);
			Vector2D Scale(double s);
			Vector2D Scale(double sX, double sY);
			Vector2D Rotate(double angle);

			double GetAngle() const;
			double GetLen() const;

			Primitives::Float2 ApplyVector(Primitives::Float2 p);
		};

		EXPORT Vector3D GetRelativeVec(const Vector3D& vec1, const Vector3D& vec2);
		EXPORT Vector2D VectorAdd(Vector2D v1, Vector2D v2);
		EXPORT Vector2D VectorSub(Vector2D v1, Vector2D v2);
		EXPORT double VectorGetAngleDifference(Vector2D v1, Vector2D v2);
		EXPORT double VectorDotProduct(Vector2D v1, Vector2D v2);
		EXPORT Vector3D VectorCrossProduct(Vector2D v1, Vector2D v2);
		EXPORT double VectorGetDeterminant(Vector2D v1, Vector2D v2);
		EXPORT void PrintProperties(Vector2D v);
		EXPORT void PrintProperties(Vector3D v);
	}

	namespace Primitives
	{
		EXPORT Float2 Transform(Float2 point, Float2 origin, Matrices::MatrixF* transform);
		EXPORT Line2D Transform(Line2D line, Float2 origin, Matrices::MatrixF* transform);
		EXPORT Float3 RotatePoint(Float3 point, const Complex::Quaternion& quat);
	}

	namespace Utility
	{
		struct EXPORT Intersect
		{
			bool isIntersecting;
			bool isCollinear;
			Primitives::Float2 pos;
		};

		EXPORT double Deg2Rad(double deg);
		EXPORT double Rad2Deg(double rad);
		EXPORT double GetSlope(Primitives::Line2D line);
		EXPORT bool IsParallel(Primitives::Line2D l1, Primitives::Line2D l2);
		EXPORT Intersect CheckIntersect(Primitives::Line2D l1, Primitives::Line2D l2);
		EXPORT Intersect CheckIntersect(Primitives::Line2D l, Primitives::Float2 p);
		EXPORT void PrintProperties(Intersect i);
		EXPORT Primitives::Line2D Vector2Line(Vectors::Vector2D vector, Primitives::Float2 pos);
		EXPORT Vectors::Vector2D Line2Vector(Primitives::Line2D line);
		EXPORT Vectors::Vector2D Line2Vector(Primitives::Float2 p1, Primitives::Float2 p2);
		EXPORT double GetDistance(Primitives::Float2 p1, Primitives::Float2 p2);
		EXPORT double GetDistance(Primitives::Float3 p1, Primitives::Float3 p2);
		EXPORT int RandInt(int min, int max);
		EXPORT inline float RoundFTo(const float& num, const int& decimal)
		{
			return roundf(pow(10, decimal) * num) / pow(10, decimal);
		}
		EXPORT Matrices::MatrixF Vec2Mat(const Vectors::Vector3D& vec, const int& front = 0);
		EXPORT Vectors::Vector3D Mat2Vec(const Matrices::MatrixF& mat, const int& front = 0);
	}
}
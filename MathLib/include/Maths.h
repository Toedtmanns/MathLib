#pragma once
#define MATHLIB

#include <math.h>
#include <stdio.h>
#include <vector>
#include <stdexcept>
#include <stdlib.h>

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
			imaginaryBase operator+(const imaginaryBase& other);
			imaginaryBase operator-(const imaginaryBase& other);
			double operator*(const imaginaryBase& other);
			double operator/(const imaginaryBase& other);
			bool operator==(const imaginaryBase& other);
			bool operator!=(const imaginaryBase& other);
			bool operator>=(const imaginaryBase& other);
			bool operator<=(const imaginaryBase& other);
			bool operator>(const imaginaryBase& other);
			bool operator<(const imaginaryBase& other);
			imaginaryBase& operator++();
			imaginaryBase& operator--();
			imaginaryBase& operator++(int);
			imaginaryBase& operator--(int);
		};

		template <typename T>
		EXPORT imaginaryBase operator*(const T& num1, const imaginaryBase& num2)
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

			double operator*(const imagI& other);
			imagK operator*(const imagJ& jNum);
			imagJ operator*(const imagK& jNum);
		};
		class EXPORT imagJ : public imaginaryBase
		{
		public:
			imagJ();
			imagJ(const double& num);
			imagJ(const imaginaryBase& base);

			double operator*(const imagJ& other);
			imagI operator*(const imagK& jNum);
			imagK operator*(const imagI& jNum);
		};
		class EXPORT imagK : public imaginaryBase
		{
		public:
			imagK();
			imagK(const double& num);
			imagK(const imaginaryBase& base);

			double operator*(const imagK& other);
			imagJ operator*(const imagI& jNum);
			imagI operator*(const imagJ& jNum);
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

			Quaternion operator+(const Quaternion& other);
			Quaternion operator*(const Quaternion& other);
			Quaternion RotateQuaternion(const Quaternion& quat);
			Primitives::Float3 RotatePoint(const Primitives::Float3& point);

			Quaternion GetInverse() const;
			Primitives::Float3 GetPoint() const;
		};

		EXPORT Quaternion QuaternionRotation(const double& angle, const double& iAxis, const double& jAxis, const double& kAxis);
		EXPORT Quaternion QuaternionRotation(const double& angle, const Primitives::Float3& axis);
		EXPORT void PrintProperties(const Quaternion& quat);
	}

	namespace Matrices
	{
		class EXPORT RowI
		{
			std::vector<int> m_Row;
			int m_Length;
		public:

			RowI(int length, int val);
			RowI(std::vector<int> row);

			void SetRow(std::vector<int>* row);
			void SetNum(int collumn, int value);

			int GetLength() const;
			int GetAt(int index) const;
			std::vector<int> GetRow() const;
		};
		class EXPORT RowF
		{
			std::vector<double> m_Row;
			int m_Length;
		public:

			RowF(int length, double val);
			RowF(std::vector<double> row);

			void SetRow(std::vector<double>* row);
			void SetNum(int collumn, double value);

			int GetLength() const;
			double GetAt(int index) const;
			std::vector<int> GetRow() const;
		};

		class EXPORT MatrixI
		{
			std::vector<RowI> m_Matrix;
			int m_Rows, m_Collumns;
		public:

			MatrixI(int rows, int collumns);
			MatrixI(std::vector<RowI> content);

			void SetRow(int index, RowI row);
			void SetNum(int row, int collumn, int value);

			int GetRowCount() const;
			int GetCollumnCount() const;
			RowI GetRow(int index) const;

			~MatrixI();
		};
		class EXPORT MatrixF
		{
			std::vector<RowF> m_Matrix;
			int m_Rows, m_Collumns;
		public:

			MatrixF(int dim);
			MatrixF(int rows, int collumns);
			MatrixF(std::vector<RowF> content);

			void SetRow(int index, RowF row);
			void SetNum(int row, int collumn, double value);

			int GetRowCount() const;
			int GetCollumnCount() const;
			RowF GetRow(int index) const;

			~MatrixF();
		};

		EXPORT void PrintContent(MatrixI* mat);
		EXPORT void PrintProperties(MatrixI* mat);
		EXPORT int MatrixGetDet(MatrixI* mat);
		EXPORT MatrixI MatrixAdd(MatrixI* mat, int value);
		EXPORT MatrixI MatrixSub(MatrixI* mat, int value);
		EXPORT MatrixI MatrixMult(MatrixI* mat, int value);
		EXPORT MatrixI MatrixDiv(MatrixI* mat, int value);
		EXPORT MatrixI MatrixAdd(MatrixI* mat1, MatrixI* mat2);
		EXPORT MatrixI MatrixSub(MatrixI* mat1, MatrixI* mat2);
		EXPORT MatrixI MatrixMult(MatrixI* mat1, MatrixI* mat2);
		EXPORT MatrixI MatrixOfMinors(MatrixI* mat);
		EXPORT MatrixI MatrixOfCofactors(MatrixI* mat);
		EXPORT MatrixI MatrixAdjugate(MatrixI* mat);
		EXPORT MatrixF MatrixInverse(MatrixI* mat);

		EXPORT void PrintContent(MatrixF* mat);
		EXPORT void PrintProperties(MatrixF* mat);
		EXPORT double MatrixGetDet(MatrixF* mat);
		EXPORT MatrixF MatrixAdd(MatrixF* mat, double value);
		EXPORT MatrixF MatrixSub(MatrixF* mat, double value);
		EXPORT MatrixF MatrixMult(MatrixF* mat, double value);
		EXPORT MatrixF MatrixDiv(MatrixF* mat, double value);
		EXPORT MatrixF MatrixAdd(MatrixF* mat1, MatrixF* mat2);
		EXPORT MatrixF MatrixSub(MatrixF* mat1, MatrixF* mat2);
		EXPORT MatrixF MatrixMult(MatrixF* mat1, MatrixF* mat2);
		EXPORT MatrixF MatrixOfMinors(MatrixF* mat);
		EXPORT MatrixF MatrixOfCofactors(MatrixF* mat);
		EXPORT MatrixF MatrixAdjugate(MatrixF* mat);
		EXPORT MatrixF MatrixInverse(MatrixF* mat);

		// Quaternion operations

		EXPORT MatrixF MatrixRotate(MatrixF* mat, const Complex::Quaternion& quat);

		// Conversion and checking functions

		EXPORT MatrixI MatrixF2I(MatrixF* mat);
		EXPORT MatrixF MatrixI2F(MatrixI* mat);
		EXPORT bool MatrixIsSquare(const MatrixI* mat, int dimension = -1);
		EXPORT bool MatrixIsSquare(const MatrixF* mat, int dimension = -1);

		// Helper functions

		EXPORT MatrixI* MatrixI2x2(int r1c1, int r1c2, int r2c1, int r2c2);
		EXPORT MatrixI* TransformI2x2();
		EXPORT MatrixI* TransformI2x2(Primitives::Float2 c1, Primitives::Float2 c2);
		EXPORT MatrixF* MatrixF2x2(double r1c1, double r1c2, double r2c1, double r2c2);
		EXPORT MatrixF* TransformF2x2();
		EXPORT MatrixF* TransformF2x2(Primitives::Float2 c1, Primitives::Float2 c2);
	}

	namespace Vectors
	{
		class EXPORT Vector3D
		{
		public:
			Primitives::Float3 direction;

			Vector3D();
			Vector3D(Primitives::Float3 dir);

			Vector3D Transform(Matrices::MatrixF* mat);
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

			Vector2D Transform(Matrices::MatrixF* mat);
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
		EXPORT Vectors::Vector3D Mat2Vec(Matrices::MatrixF* mat, const int& front = 0);
	}
}
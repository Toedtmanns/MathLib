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
	namespace Complex
	{
		struct imaginaryNum
		{
			double num;

			imaginaryNum();
			imaginaryNum(double num);

			template <typename T>
			imaginaryNum operator*(const T &other)
			{
				return imaginaryNum(num * other);
			}
			double operator+(const imaginaryNum &other);
			double operator-(const imaginaryNum &other);
			double operator*(const imaginaryNum &other);
			double operator/(const imaginaryNum &other);
			bool operator==(const imaginaryNum &other);
			bool operator!=(const imaginaryNum &other);
			bool operator>=(const imaginaryNum &other);
			bool operator<=(const imaginaryNum &other);
			bool operator>(const imaginaryNum &other);
			bool operator<(const imaginaryNum &other);
			imaginaryNum& operator++();
			imaginaryNum& operator--();
			double operator++(int);
			double operator--(int);
		};

		template <typename T>
		imaginaryNum operator*(const T &num1, const imaginaryNum &num2)
		{
			return imaginaryNum(num1 * num2.num);
		}

		class Quaternion
		{
			double real;
			imaginaryNum i, j, k;

			Quaternion();
			Quaternion(double real, imaginaryNum i, imaginaryNum j, imaginaryNum k);

			Quaternion operator+(const Quaternion &other);
			Quaternion operator*(const Quaternion &other);
		};
	}

	typedef Complex::imaginaryNum imaginary;

	namespace Primitives
	{
		struct EXPORT Point3D
		{
			double x;
			double y;
			double z;

			Point3D()
			{
				x = 0;
				y = 0;
				z = 0;
			}
			Point3D(double x, double y, double z)
			{
				this->x = x;
				this->y = y;
				this->z = z;
			}
		};

		struct EXPORT Point2D
		{
			double x;
			double y;

			Point2D()
			{
				x = 0;
				y = 0;
			}
			Point2D(double x, double y)
			{
				this->x = x;
				this->y = y;
			}
		};

		struct EXPORT Line2D
		{
			Point2D p1;
			Point2D p2;
			float normal;

			Line2D();
			Line2D(double x1, double y1, double x2, double y2);
			Line2D(Point2D p1, Point2D p2);
		};

		EXPORT bool operator==(const Point3D &p1, const Point3D &p2);
		EXPORT bool operator!=(const Point3D &p1, const Point3D &p2);
		EXPORT Point2D PointAdd(const Point2D &p1, const Point2D &p2);
		EXPORT Point2D operator+(const Point2D &p1, const Point2D &p2);
		EXPORT Point3D PointAdd(const Point3D &p1, const Point3D &p2);
		EXPORT Point3D operator+(const Point3D &p1, const Point3D &p2);
		EXPORT void PrintProperties(Point2D p);
		EXPORT void PrintProperties(Line2D l);
	}

	namespace Matrices
	{
		class EXPORT RowI
		{
			std::vector<int> row;
			int length;
		public:

			RowI(int length, int val);
			RowI(std::vector<int> row);

			int GetLength() const;
			int GetAt(int index);
			std::vector<int> GetRow() const;
			void SetRow(std::vector<int> *row);
			void SetNum(int collumn, int value);
		};
		class EXPORT RowF
		{
			std::vector<double> row;
			int length;
		public:

			RowF(int length, double val);
			RowF(std::vector<double> row);

			int GetLength() const;
			double GetAt(int index);
			std::vector<int> GetRow() const;
			void SetRow(std::vector<double> *row);
			void SetNum(int collumn, double value);
		};

		class EXPORT MatrixI
		{
			std::vector<RowI> matrix;
			int rows, collumns;
		public:

			MatrixI(int rows, int collumns);
			MatrixI(std::vector<RowI> content);

			int GetRowCount() const;
			int GetCollumnCount() const;
			RowI GetRow(int index);
			void SetRow(int index, RowI row);
			void SetNum(int row, int collumn, int value);
		};
		class EXPORT MatrixF
		{
			std::vector<RowF> matrix;
			int rows, collumns;
		public:

			MatrixF(int dim);
			MatrixF(int rows, int collumns);
			MatrixF(std::vector<RowF> content);

			int GetRowCount() const;
			int GetCollumnCount() const;
			RowF GetRow(int index);
			void SetRow(int index, RowF row);
			void SetNum(int row, int collumn, double value);

			~MatrixF();
		};

		EXPORT void PrintContent(MatrixI *mat);
		EXPORT void PrintProperties(MatrixI *mat);
		EXPORT int MatrixGetDet(MatrixI *mat);
		EXPORT MatrixI MatrixAdd(MatrixI *mat, int value);
		EXPORT MatrixI MatrixSub(MatrixI *mat, int value);
		EXPORT MatrixI MatrixMult(MatrixI *mat, int value);
		EXPORT MatrixI MatrixDiv(MatrixI *mat, int value);
		EXPORT MatrixI MatrixAdd(MatrixI *mat1, MatrixI *mat2);
		EXPORT MatrixI MatrixSub(MatrixI *mat1, MatrixI *mat2);
		EXPORT MatrixI MatrixMult(MatrixI *mat1, MatrixI *mat2);
		EXPORT MatrixI MatrixOfMinors(MatrixI *mat);
		EXPORT MatrixI MatrixOfCofactors(MatrixI *mat);
		EXPORT MatrixI MatrixAdjugate(MatrixI *mat);
		EXPORT MatrixF MatrixInverse(MatrixI *mat);

		EXPORT void PrintContent(MatrixF *mat);
		EXPORT void PrintProperties(MatrixF *mat);
		EXPORT double MatrixGetDet(MatrixF *mat);
		EXPORT MatrixF MatrixAdd(MatrixF *mat, double value);
		EXPORT MatrixF MatrixSub(MatrixF *mat, double value);
		EXPORT MatrixF MatrixMult(MatrixF *mat, double value);
		EXPORT MatrixF MatrixDiv(MatrixF *mat, double value);
		EXPORT MatrixF MatrixAdd(MatrixF *mat1, MatrixF *mat2);
		EXPORT MatrixF MatrixSub(MatrixF *mat1, MatrixF *mat2);
		EXPORT MatrixF MatrixMult(MatrixF *mat1, MatrixF *mat2);
		EXPORT MatrixF MatrixOfMinors(MatrixF *mat);
		EXPORT MatrixF MatrixOfCofactors(MatrixF *mat);
		EXPORT MatrixF MatrixAdjugate(MatrixF *mat);
		EXPORT MatrixF MatrixInverse(MatrixF *mat);

		EXPORT MatrixI MatrixF2I(MatrixF *mat);
		EXPORT MatrixF MatrixI2F(MatrixI *mat);
		EXPORT bool MatrixIsSquare(const MatrixI &mat, int dimension = -1);
		EXPORT bool MatrixIsSquare(const MatrixF &mat, int dimension = -1);

		// Helper functions

		EXPORT MatrixI *MatrixI2x2(int r1c1, int r1c2, int r2c1, int r2c2);
		EXPORT MatrixI *TransformI2x2();
		EXPORT MatrixI *TransformI2x2(Primitives::Point2D c1, Primitives::Point2D c2);
		EXPORT MatrixF *MatrixF2x2(double r1c1, double r1c2, double r2c1, double r2c2);
		EXPORT MatrixF *TransformF2x2();
		EXPORT MatrixF *TransformF2x2(Primitives::Point2D c1, Primitives::Point2D c2);
	}

	namespace Vectors
	{
		class EXPORT Vector3D
		{
		public:
			Primitives::Point3D direction;

			Vector3D();
			Vector3D(Primitives::Point3D dir);

			Vector3D Transform(Matrices::MatrixF *mat);
			Vector3D Scale(double s);
			Vector3D Scale(double sX, double sY, double sZ);
			Vector3D Rotate(double angle, int axis);
			Vector3D Rotate(double angle, Vector3D *normal);

			double GetAngle(int axis);
			double GetAngle(Vector3D *normal);
			double GetLen();
		};

		/// <summary>
		/// A class which handles vectors and vector calculations
		/// </summary>
		class EXPORT Vector2D
		{
		public:
			Primitives::Point2D direction;

			Vector2D();
			Vector2D(Primitives::Point2D dir);
			Vector2D(double angle, double length);

			Vector2D Transform(Matrices::MatrixF *mat);
			Vector2D Scale(double s);
			Vector2D Scale(double sX, double sY);
			Vector2D Rotate(double angle);

			double GetAngle();
			double GetLen();

			Primitives::Point2D ApplyVector(Primitives::Point2D p);
		};

		EXPORT Vector3D GetRelativeVec(const Vector3D &vec1, const Vector3D &vec2);
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
		EXPORT Point2D Transform(Point2D point, Point2D origin, Matrices::MatrixF *transform);
		EXPORT Line2D Transform(Line2D line, Point2D origin, Matrices::MatrixF *transform);
	}

	namespace Utility
	{
		struct EXPORT Intersect
		{
			bool isIntersecting;
			bool isCollinear;
			Primitives::Point2D pos;
		};

		EXPORT double Deg2Rad(double deg);
		EXPORT double Rad2Deg(double rad);
		EXPORT double GetSlope(Primitives::Line2D line);
		EXPORT bool IsParallel(Primitives::Line2D l1, Primitives::Line2D l2);
		EXPORT Intersect CheckIntersect(Primitives::Line2D l1, Primitives::Line2D l2);
		EXPORT Intersect CheckIntersect(Primitives::Line2D l, Primitives::Point2D p);
		EXPORT void PrintProperties(Intersect i);
		EXPORT Primitives::Line2D Vector2Line(Vectors::Vector2D vector, Primitives::Point2D pos);
		EXPORT Vectors::Vector2D Line2Vector(Primitives::Line2D line);
		EXPORT Vectors::Vector2D Line2Vector(Primitives::Point2D p1, Primitives::Point2D p2);
		EXPORT double GetDistance(Primitives::Point2D p1, Primitives::Point2D p2);
		EXPORT double GetDistance(Primitives::Point3D p1, Primitives::Point3D p2);
		EXPORT int RandInt(int min, int max);
		EXPORT inline float RoundFTo(const float &num, const int &decimal)
		{
			return roundf(pow(10, decimal) * num) / pow(10, decimal);
		}
		EXPORT Matrices::MatrixF Vec2Mat(const Vectors::Vector3D &vec, const int &front = 0);
		EXPORT Vectors::Vector3D Mat2Vec(Matrices::MatrixF *mat, const int &front = 0);
	}
}
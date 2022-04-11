#include "../include/Maths.h"
#include <stdio.h>
#include <cmath>

namespace MathLib
{
	Float4::Float4()
			: x(0), y(0), z(0), w(0)
		{

		}
	Float4::Float4(const double& x, const double& y, const double& z, const double& w)
			: x(x), y(y), z(z), w(w)
		{

		}
	Float4::Float4(const Int4& other)
		{
			x = other.x;
			y = other.y;
			z = other.z;
			w = other.w;
		}
	void Float4::operator=(const Int4& other)
		{
			x = other.x;
			y = other.y;
			z = other.z;
			w = other.w;
		}
	double& Float4::operator[](const unsigned int& index)
		{
			return content[index];
		}
	const double& Float4::operator[](const unsigned int& index) const
		{
			return content[index];
		}

	Float3::Float3()
			: x(0), y(0), z(0)
		{

		}
	Float3::Float3(const double& x, const double& y, const double& z)
			: x(x), y(y), z(z)
		{

		}
	Float3::Float3(const Int3& other)
		{
			x = other.x;
			y = other.y;
			z = other.z;
		}
	void Float3::operator=(const Int3& other)
		{
			x = other.x;
			y = other.y;
			z = other.z;
		}
	double& Float3::operator[](const unsigned int& index)
		{
			
			return content[index];
			
		}
	const double& Float3::operator[](const unsigned int& index) const
		{
			return content[index];
		}

	Float2::Float2()
			: x(0), y(0)
		{

		}
	Float2::Float2(const double& x, const double& y)
			: x(x), y(y)
		{

		}
	Float2::Float2(const Int2& other)
		{
			x = other.x;
			y = other.y;
		}
	void Float2::operator=(const Int2& other)
		{
			x = other.x;
			y = other.y;
		}
	double& Float2::operator[](const unsigned int& index)
		{
			return content[index];
		}
	const double& Float2::operator[](const unsigned int& index) const
		{
			return content[index];
		}

	Int4::Int4()
			: x(0), y(0), z(0), w(0)
		{

		}
	Int4::Int4(const int& x, const int& y, const int& z, const int& w)
			: x(x), y(y), z(z), w(w)
		{

		}
	Int4::Int4(const Float4& other)
		{
			x = (int) round(other.x);
			y = (int) round(other.y);
			z = (int) round(other.z);
			w = (int) round(other.w);
		}
	void Int4::operator=(const Float4& other)
		{
			x = (int) round(other.x);
			y = (int) round(other.y);
			z = (int) round(other.z);
			w = (int) round(other.w);
		}
	int& Int4::operator[](const unsigned int& index)
		{
			return content[index];
		}
	const int& Int4::operator[](const unsigned int& index) const
		{
			return content[index];
		}

	Int3::Int3()
			: x(0), y(0), z(0)
		{

		}
	Int3::Int3(const int& x, const int& y, const int& z)
			: x(x), y(y), z(z)
		{

		}
	Int3::Int3(const Float3& other)
		{
			x = (int) round(other.x);
			y = (int) round(other.y);
			z = (int) round(other.z);
		}
	void Int3::operator=(const Float3& other)
		{
			x = (int) round(other.x);
			y = (int) round(other.y);
			z = (int) round(other.z);
		}
	int& Int3::operator[](const unsigned int& index)
		{
			
			return content[index];
			
		}
	const int& Int3::operator[](const unsigned int& index) const
		{
			return content[index];
		}

	Int2::Int2()
			: x(0), y(0)
		{

		}
	Int2::Int2(const int& x, const int& y)
			: x(x), y(y)
		{

		}
	Int2::Int2(const Float2& other)
		{
			x = (int) round(other.x);
			y = (int) round(other.y);
		}
	void Int2::operator=(const Float2& other)
		{
			x = (int) round(other.x);
			y = (int) round(other.y);
		}
	int& Int2::operator[](const unsigned int& index)
		{
			return content[index];
		}
	const int& Int2::operator[](const unsigned int& index) const
		{
			return content[index];
		}

	bool operator==(const Float2& f1, const Float2& f2)
		{
			return f1.x == f2.x && f1.y == f2.y;
		}
	bool operator!=(const Float2& f1, const Float2& f2)
		{
			return f1.x != f2.x || f1.y != f2.y;
		}
	bool operator==(const Float3& f1, const Float3& f2)
		{
			return f1.x == f2.x && f1.y == f2.y && f1.z == f2.z;
		}
	bool operator!=(const Float3& f1, const Float3& f2)
		{
			return f1.x != f2.x || f1.y != f2.y || f1.z != f2.z;
		}
	bool operator==(const Float4& f1, const Float4& f2)
		{
			return f1.x == f2.x && f1.y == f2.y && f1.z == f2.z && f1.w == f2.w;
		}
	bool operator!=(const Float4& f1, const Float4& f2)
		{
			return f1.x != f2.x || f1.y != f2.y || f1.z != f2.z || f1.w != f2.w;
		}
	Float2 operator+(const Float2& f1, const Float2& f2)
		{
			return Float2(f1.x + f2.x, f1.y + f2.y);
		}
	Float3 operator+(const Float3& f1, const Float3& f2)
		{
			return Float3(f1.x + f2.x, f1.y + f2.y, f1.z + f2.z);
		}
	Float4 operator+(const Float4& f1, const Float4& f2)
		{
			return Float4(f1.x + f2.x, f1.y + f2.y, f1.z + f2.z, f1.w + f2.w);
		}
	Float2 operator-(const Float2& f1, const Float2& f2)
		{
			return Float2(f1.x - f2.x, f1.y - f2.y);
		}
	Float3 operator-(const Float3& f1, const Float3& f2)
		{
			return Float3(f1.x - f2.x, f1.y - f2.y, f1.z - f2.z);
		}
	Float4 operator-(const Float4& f1, const Float4& f2)
		{
			return Float4(f1.x - f2.x, f1.y - f2.y, f1.z - f2.z, f1.w - f2.w);
		}
	void operator+=(Float2& f1, const Float2& f2)
		{
			f1 = Float2(f1.x + f2.x, f1.y + f2.y);
		}
	void operator+=(Float3& f1, const Float3& f2)
		{
			f1 = Float3(f1.x + f2.x, f1.y + f2.y, f1.z + f2.z);
		}
	void operator+=(Float4& f1, const Float4& f2)
		{
			f1 = Float4(f1.x + f2.x, f1.y + f2.y, f1.z + f2.z, f1.w + f2.w);
		}
	void operator-=(Float2& f1, const Float2& f2)
		{
			f1 = Float2(f1.x - f2.x, f1.y - f2.y);
		}
	void operator-=(Float3& f1, const Float3& f2)
		{
			f1 = Float3(f1.x - f2.x, f1.y - f2.y, f1.z - f2.z);
		}
	void operator-=(Float4& f1, const Float4& f2)
		{
			f1 = Float4(f1.x - f2.x, f1.y - f2.y, f1.z - f2.z, f1.w - f2.w);
		}
	Float2 operator*(const Float2& f1, const Float2& f2)
		{
			return Float2(f1.x * f2.x, f1.y * f2.y);
		}
	Float3 operator*(const Float3& f1, const Float3& f2)
		{
			return Float3(f1.x * f2.x, f1.y * f2.y, f1.z * f2.z);
		}
	Float4 operator*(const Float4& f1, const Float4& f2)
		{
			return Float4(f1.x * f2.x, f1.y * f2.y, f1.z * f2.z, f1.w * f2.w);
		}
	Float2 operator*(const Float2& f1, const double& f2)
		{
			return Float2(f1.x * f2, f1.y * f2);
		}
	Float3 operator*(const Float3& f1, const double& f2)
		{
			return Float3(f1.x * f2, f1.y * f2, f1.z * f2);
		}
	Float4 operator*(const Float4& f1, const double& f2)
		{
			return Float4(f1.x * f2, f1.y * f2, f1.z * f2, f1.w * f2);
		}

	bool operator==(const Int2& i1, const Int2& i2)
		{
			return i1.x == i2.x && i1.y == i2.y;
		}
	bool operator!=(const Int2& i1, const Int2& i2)
		{
			return i1.x != i2.x || i1.y != i2.y;
		}
	bool operator==(const Int3& i1, const Int3& i2)
		{
			return i1.x == i2.x && i1.y == i2.y && i1.z == i2.z;
		}
	bool operator!=(const Int3& i1, const Int3& i2)
		{
			return i1.x != i2.x || i1.y != i2.y || i1.z != i2.z;
		}
	bool operator==(const Int4& i1, const Int4& i2)
		{
			return i1.x == i2.x && i1.y == i2.y && i1.z == i2.z && i1.w == i2.w;
		}
	bool operator!=(const Int4& i1, const Int4& i2)
		{
			return i1.x != i2.x || i1.y != i2.y || i1.z != i2.z || i1.w != i2.w;
		}
	Int2 operator+(const Int2& i1, const Int2& i2)
		{
			return Int2(i1.x + i2.x, i1.y + i2.y);
		}
	Int3 operator+(const Int3& i1, const Int3& i2)
		{
			return Int3(i1.x + i2.x, i1.y + i2.y, i1.z + i2.z);
		}
	Int4 operator+(const Int4& i1, const Int4& i2)
		{
			return Int4(i1.x + i2.x, i1.y + i2.y, i1.z + i2.z, i1.w + i2.w);
		}
	Int2 operator-(const Int2& i1, const Int2& i2)
		{
			return Int2(i1.x - i2.x, i1.y - i2.y);
		}
	Int3 operator-(const Int3& i1, const Int3& i2)
		{
			return Int3(i1.x - i2.x, i1.y - i2.y, i1.z - i2.z);
		}
	Int4 operator-(const Int4& i1, const Int4& i2)
		{
			return Int4(i1.x - i2.x, i1.y - i2.y, i1.z - i2.z, i1.w - i2.w);
		}
	void operator+=(Int2& i1, const Int2& i2)
		{
			i1 = Int2(i1.x + i2.x, i1.y + i2.y);
		}
	void operator+=(Int3& i1, const Int3& i2)
		{
			i1 = Int3(i1.x + i2.x, i1.y + i2.y, i1.z + i2.z);
		}
	void operator+=(Int4& i1, const Int4& i2)
		{
			i1 = Int4(i1.x + i2.x, i1.y + i2.y, i1.z + i2.z, i1.w + i2.w);
		}
	void operator-=(Int2& i1, const Int2& i2)
		{
			i1 = Int2(i1.x - i2.x, i1.y - i2.y);
		}
	void operator-=(Int3& i1, const Int3& i2)
		{
			i1 = Int3(i1.x - i2.x, i1.y - i2.y, i1.z - i2.z);
		}
	void operator-=(Int4& i1, const Int4& i2)
		{
			i1 = Int4(i1.x - i2.x, i1.y - i2.y, i1.z - i2.z, i1.w - i2.w);
		}
	Int2 operator*(const Int2& i1, const Int2& i2)
		{
			return Int2(i1.x * i2.x, i1.y * i2.y);
		}
	Int3 operator*(const Int3& i1, const Int3& i2)
		{
			return Int3(i1.x * i2.x, i1.y * i2.y, i1.z * i2.z);
		}
	Int4 operator*(const Int4& i1, const Int4& i2)
		{
			return Int4(i1.x * i2.x, i1.y * i2.y, i1.z * i2.z, i1.w * i2.w);
		}
	Int2 operator*(const Int2& i1, const int& i2)
		{
			return Int2(i1.x * i2, i1.y * i2);
		}
	Int3 operator*(const Int3& i1, const int& i2)
		{
			return Int3(i1.x * i2, i1.y * i2, i1.z * i2);
		}
	Int4 operator*(const Int4& i1, const int& i2)
		{
			return Int4(i1.x * i2, i1.y * i2, i1.z * i2, i1.w * i2);
		}

	Float2 Transform(Float2 point, Float2 origin, const Matrices::MatrixF& transform)
		{
			Vectors::Vector2D vec = Utility::Line2Vector(origin, point);
			vec.Transform(transform);

			return origin + vec.direction;
		}
	Line2D Transform(Line2D line, Float2 origin, const Matrices::MatrixF& transform)
		{
			Vectors::Vector2D vec1 = Utility::Line2Vector(origin, line.p1);
			Vectors::Vector2D vec2 = Utility::Line2Vector(origin, line.p2);
			vec1.Transform(transform);
			vec2.Transform(transform);

			return Line2D(vec1.direction + origin, vec2.direction + origin);
		}
	Float3 RotatePoint(const Float3& point, Complex::Quaternion& quat)
		{
			Complex::Quaternion retQuat(point);
			retQuat = quat * retQuat * quat.GetInverse();
			return retQuat.GetPoint();
		}

	void PrintProperties(const Float2& p)
		{
			printf("X: %f, Y: %f\n", p.x, p.y);
		}
	void PrintProperties(const Float3& p)
		{
			printf("X: %.3f, Y: %.3f, Z: %.3f\n", p.x, p.y, p.z);
		}

	// Lines

	Line2D::Line2D()
		{
			p1 = Float2();
			p2 = Float2();
			normal = 0;
		}
	Line2D::Line2D(const double& x1, const double& y1, const double& x2, const double& y2)
		{
			p1.x = x1;
			p1.y = y1;
			p2.x = x2;
			p2.y = y2;
			normal = Utility::Line2Vector(Float2(x1, y1), Float2(x2, y2)).GetAngle() + 90;
		}
	Line2D::Line2D(const Float2& p1, const Float2& p2)
		{
			this->p1.x = p1.x;
			this->p1.y = p1.y;
			this->p2.x = p2.x;
			this->p2.y = p2.y;

			normal = Utility::Line2Vector(p1, p2).GetAngle() + 90;
			if (normal < 0)
				normal = 360 + normal;
		}
	bool Line2D::operator==(const Line2D& other)
		{
			return (p1 == other.p1 && p2 == other.p2);
		}
	bool Line2D::operator!=(const Line2D& other)
		{
			return !(p1 == other.p1 && p2 == other.p2);
		}

	// Intersects

	Intersect::Intersect()
			: pos({nan(""), nan("")}), isIntersecting(false), isCollinear(false), line1(), line2()
		{

		}
	Intersect::Intersect(const Float2& pos, const bool& intersecting, const bool& parallel)
			: pos(pos), isIntersecting(intersecting), isCollinear(parallel), line1(), line2()
		{

		}

	// Other functions

	double GetLength(const Line2D& line)
		{
			return sqrt(pow(line.p2.x - line.p1.x, 2) * pow(line.p2.y - line.p1.y, 2));
		}
	double GetDistance(const Float2& p1, const Float2& p2)
		{
			return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
		}
	double GetDistance(const Float3& p1, const Float3& p2)
		{
			return sqrt(pow(GetDistance({p1.x, p1.y}, {p2.x, p2.y}), 2) + pow(p2.z - p1.z, 2));
		}

	double GetSlope(const Line2D& line)
		{
			return (line.p2.y - line.p1.y) / (line.p2.x - line.p1.x);
		}
	double GetYAxisSection(const Line2D& line)
		{
			return line.p1.y - GetSlope(line) * line.p1.x;
		}
	bool IsParallel(const Line2D& l1, const Line2D& l2)
		{
			double slopeL1 = GetSlope(l1);
			double slopeL2 = GetSlope(l2);

			if (slopeL1 == slopeL2 || isinf(slopeL1) && isinf(slopeL2))
				return true;
			else
				return false;
		}
	Float2 GetIntersectPos(const Line2D& l1, const Line2D& l2)
		{
			Float2 iPos(
				(GetYAxisSection(l2) - GetYAxisSection(l1)) / (GetSlope(l1) - GetSlope(l2)),
				GetSlope(l1) * ((GetYAxisSection(l2) - GetYAxisSection(l1)) / (GetSlope(l1) - GetSlope(l2))) + GetYAxisSection(l1)
			);

			return iPos;
		}
	Intersect GetIntersect(const Line2D& l1, const Line2D& l2)
		{
			double m1 = GetSlope(l1);
			double m2 = GetSlope(l2);
			double b1 = l1.p1.y - l1.p1.x * m1;
			double b2 = l2.p2.y - l2.p2.x * m2;
			Intersect intersect = Intersect();

			intersect.pos.x = (b2 - b1) / (m1 - m2) * 100 / 100;
			intersect.pos.y = (m1 * intersect.pos.x + b1) * 100 / 100;
			intersect.line1 = l1;
			intersect.line2 = l2;

			if (isnan(intersect.pos.x))
			{
				intersect.isCollinear = true;

				// Check intersect if one line goes straight up or down and find the intersect

				if (isinf(m1) && l1.p1.x <= fmax(l2.p1.x, l2.p2.x) && l1.p1.x >= fmin(l2.p1.x, l2.p2.x))
				{
					Float2 possInt = Float2(l1.p1.x, m2 * l1.p1.x + b2);
					if (GetIntersect(l1, possInt).isIntersecting)
					{
						intersect.pos = possInt;
						intersect.isIntersecting = true;
						intersect.isCollinear = false;
					}
				}
				else if (isinf(m2) && l2.p1.x <= fmax(l1.p1.x, l1.p2.x) && l2.p1.x >= fmin(l1.p1.x, l1.p2.x))
				{
					Float2 possInt = Float2(l2.p1.x, m1 * l2.p1.x + b1);
					if (GetIntersect(l2, possInt).isIntersecting)
					{
						intersect.pos = possInt;
						intersect.isIntersecting = true;
						intersect.isCollinear = false;
					}
				}

				// Check intersect if both lines go straight up or down

				else if (GetIntersect(l1, l2.p1).isIntersecting)
					intersect.isIntersecting = true;
				else if (GetIntersect(l1, l2.p2).isIntersecting)
					intersect.isIntersecting = true;
				else if (GetIntersect(l2, l1.p1).isIntersecting)
					intersect.isIntersecting = true;
				else if (GetIntersect(l2, l1.p2).isIntersecting)
					intersect.isIntersecting = true;
			}
			else if (isfinite(intersect.pos.x))
			{
				if ((intersect.pos.x >= fmin(l1.p1.x, l1.p2.x) && intersect.pos.x <= fmax(l1.p1.x, l1.p2.x)) && (intersect.pos.x >= fmin(l2.p1.x, l2.p2.x) && intersect.pos.x <= fmax(l2.p1.x, l2.p2.x)))
					intersect.isIntersecting = true;
			}

			return intersect;
		}
	Intersect GetIntersect(const Line2D& l, const Float2& p)
		{
			double m = GetSlope(l);
			double b = l.p1.y - l.p1.x * m;
			Intersect intersect = Intersect();

			if (m * p.x + b == p.y)
			{
				intersect.isCollinear = true;
				if (l.p1.x <= p.x && p.x <= l.p2.x)
				{
					intersect.isIntersecting = true;
					intersect.pos.x = p.x;
					intersect.pos.y = m * p.x + b;
				}
			}
			else if (isinf(m) && p.x == l.p1.x)
			{
				intersect.isCollinear = true;
				if (p.y >= fmin(l.p1.y, l.p2.y) && p.y <= fmax(l.p1.y, l.p2.y))
				{
					intersect.isIntersecting = true;
					intersect.pos = {p.x, p.y};
				}
			}

			return intersect;
		}

	// Matrix calculations

	Float2 operator*(const Float2& point, const Matrices::MatrixF& matrix)
		{
			Float2 retFloat = Float2(
				point.x * matrix.GetNum(0, 0) + point.y * matrix.GetNum(1, 0),
				point.x * matrix.GetNum(0, 1) + point.y * matrix.GetNum(1, 1));
			return retFloat;
		}

	void PrintProperties(const Line2D& l)
		{
			printf("X1: %f, Y1: %f, X2: %f, Y2: %f\n", l.p1.x, l.p1.y, l.p2.x, l.p2.y);
		}
	void PrintProperties(const Intersect& i)
		{
			printf("X: %.3f, Y: %.3f, Intersecting: %d, Collinear: %d\n", i.pos.x, i.pos.y, i.isIntersecting, i.isCollinear);
		}
}
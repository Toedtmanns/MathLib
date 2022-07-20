#include "../include/MathLib.hpp"
#include <stdio.h>
#include <math.h>

namespace MathLib
{
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

	void PrintProperties(const Float2& p)
	{
		printf("X: %.3f, Y: %.3f\n", p.x, p.y);
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
		normal = Vector2D(Float2(x1, y1), Float2(x2, y2)).GetAngle() + 90;
	}
	Line2D::Line2D(const Float2& p1, const Float2& p2)
	{
		this->p1.x = p1.x;
		this->p1.y = p1.y;
		this->p2.x = p2.x;
		this->p2.y = p2.y;

		normal = Vector2D(p1, p2).GetAngle() + 90;
		if (normal < 0)
			normal += 360;
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
		: pos(nan(""), nan("")), isIntersecting(false), isCollinear(false), line1(), line2()
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
				intersect.pos = Float2(p.x, p.y);
			}
		}

		return intersect;
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
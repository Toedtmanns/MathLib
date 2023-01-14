#include "../include/MathLib/MathLib.hpp"
#include <stdio.h>
#include <math.h>

namespace MathLib
{
	// Intersects

	Intersect::Intersect()
		: pos(nan(""), nan("")), isIntersecting(false), isCollinear(false), line1(), line2()
	{

	}
	Intersect::Intersect(const Float2& pos, const bool intersecting, const bool parallel)
		: pos(pos), isIntersecting(intersecting), isCollinear(parallel), line1(), line2()
	{

	}

	// Other functions

	float GetSlope(const Line2D& line)
	{
		return (line.p2.y - line.p1.y) / (line.p2.x - line.p1.x);
	}
	float GetYAxisSection(const Line2D& line)
	{
		return line.p1.y - GetSlope(line) * line.p1.x;
	}
	bool IsParallel(const Line2D& l1, const Line2D& l2)
	{
		float slopeL1 = GetSlope(l1);
		float slopeL2 = GetSlope(l2);

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
		float m1 = GetSlope(l1);
		float m2 = GetSlope(l2);
		float b1 = l1.p1.y - l1.p1.x * m1;
		float b2 = l2.p2.y - l2.p2.x * m2;
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
		float m = GetSlope(l);
		float b = l.p1.y - l.p1.x * m;
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
	void PrintProperties(const Float2& p)
	{
		printf("X: %.3f, Y: %.3f\n", p.x, p.y);
	}
	void PrintProperties(const Float3& p)
	{
		printf("X: %.3f, Y: %.3f, Z: %.3f\n", p.x, p.y, p.z);
	}
}
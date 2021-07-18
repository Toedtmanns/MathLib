#include "../include/Maths.h"

namespace MathLib
{
	namespace Utility
	{
		double Deg2Rad(double deg)
		{
			return deg * PI / 180;
		}
		double Rad2Deg(double rad)
		{
			return rad * 180 / PI;
		}
		double GetSlope(Primitives::Line2D line)
		{
			return (line.p2.y - line.p1.y) / (line.p2.x - line.p1.x);
		}
		bool IsParallel(Primitives::Line2D l1, Primitives::Line2D l2)
		{
			double slopeL1 = GetSlope(l1);
			double slopeL2 = GetSlope(l2);

			if (slopeL1 == slopeL2 || isinf(slopeL1) && isinf(slopeL2))
				return true;
			else
				return false;
		}
		Intersect CheckIntersect(Primitives::Line2D l1, Primitives::Line2D l2)
		{
			double m1 = GetSlope(l1);
			double m2 = GetSlope(l2);
			double b1 = l1.p1.y - l1.p1.x * m1;
			double b2 = l2.p2.y - l2.p2.x * m2;
			Intersect intersect = {false, false, Primitives::Point2D()};

			intersect.pos.x = (b2 - b1) / (m1 - m2) * 100 / 100;
			intersect.pos.y = (m1 * intersect.pos.x + b1) * 100 / 100;

			if (isnan(intersect.pos.x))
			{
				intersect.isCollinear = true;

				// Check intersect if one line goes straight up or down and find the intersect

				if (isinf(m1) && l1.p1.x <= fmax(l2.p1.x, l2.p2.x) && l1.p1.x >= fmin(l2.p1.x, l2.p2.x))
				{
					Primitives::Point2D possInt = Primitives::Point2D(l1.p1.x, m2 * l1.p1.x + b2);
					if (CheckIntersect(l1, possInt).isIntersecting)
					{
						intersect.pos = possInt;
						intersect.isIntersecting = true;
						intersect.isCollinear = false;
					}
				}
				else if (isinf(m2) && l2.p1.x <= fmax(l1.p1.x, l1.p2.x) && l1.p1.x >= fmin(l1.p1.x, l1.p2.x))
				{
					Primitives::Point2D possInt = Primitives::Point2D(l2.p1.x, m1 * l2.p1.x + b1);
					if (CheckIntersect(l2, possInt).isIntersecting)
					{
						intersect.pos = possInt;
						intersect.isIntersecting = true;
						intersect.isCollinear = false;
					}
				}

				// Check intersect if both lines go straight up or down

				else if (CheckIntersect(l1, l2.p1).isIntersecting)
					intersect.isIntersecting = true;
				else if (CheckIntersect(l1, l2.p2).isIntersecting)
					intersect.isIntersecting = true;
				else if (CheckIntersect(l2, l1.p1).isIntersecting)
					intersect.isIntersecting = true;
				else if (CheckIntersect(l2, l1.p2).isIntersecting)
					intersect.isIntersecting = true;
			}
			else if (isfinite(intersect.pos.x))
			{
				if ((intersect.pos.x >= fmin(l1.p1.x, l1.p2.x) && intersect.pos.x <= fmax(l1.p1.x, l1.p2.x)) && (intersect.pos.x >= fmin(l2.p1.x, l2.p2.x) && intersect.pos.x <= fmax(l2.p1.x, l2.p2.x)))
					intersect.isIntersecting = true;
			}

			return intersect;
		}
		Intersect CheckIntersect(Primitives::Line2D l, Primitives::Point2D p)
		{
			double m = GetSlope(l);
			double b = l.p1.y - l.p1.x * m;
			Intersect intersect = {false, false, Primitives::Point2D()};

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
		void PrintProperties(Intersect i)
		{
			printf("X: %.3f, Y: %.3f, Intersecting: %d, Collinear: %d\n", i.pos.x, i.pos.y, i.isIntersecting, i.isCollinear);
		}
		Primitives::Line2D Vector2Line(Vectors::Vector2D vector, Primitives::Point2D pos)
		{
			Primitives::Line2D line = Primitives::Line2D();
			line.p1 = pos;
			line.p2.x = vector.direction.x + pos.x;
			line.p2.y = vector.direction.y + pos.y;

			return line;
		}
		Vectors::Vector2D Line2Vector(Primitives::Line2D line)
		{
			Vectors::Vector2D vector = Vectors::Vector2D();
			vector.direction = Primitives::Point2D(line.p2.x - line.p1.x, line.p2.y - line.p1.y);
			return vector;
		}
		Vectors::Vector2D Line2Vector(Primitives::Point2D p1, Primitives::Point2D p2)
		{
			Vectors::Vector2D vector = Vectors::Vector2D();
			vector.direction = Primitives::Point2D(p2.x - p1.x, p2.y - p1.y);
			return vector;
		}
		double GetDistance(Primitives::Point2D p1, Primitives::Point2D p2)
		{
			return sqrt(powf(p2.x - p1.x, 2) + powf(p2.y - p1.y, 2));
		}
		double GetDistance(Primitives::Point3D p1, Primitives::Point3D p2)
		{
			return sqrt(powf(GetDistance({p1.x, p1.y}, {p2.x, p2.y}), 2) + powf(p2.z - p1.z, 2));
		}
	}
}
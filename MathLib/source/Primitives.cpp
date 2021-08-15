#include "../include/Maths.h"

namespace MathLib
{
	namespace Primitives
	{
		Line2D::Line2D()
		{
			p1 = Point2D();
			p2 = Point2D();
			normal = 0;
		}			
		Line2D::Line2D(double x1, double y1, double x2, double y2)
		{
			p1.x = x1;
			p1.y = y1;
			p2.x = x2;
			p2.y = y2;
			normal = Utility::Line2Vector(Point2D(x1, y1), Point2D(x2, y2)).GetAngle() - 90;
		}
		Line2D::Line2D(Point2D p1, Point2D p2)
		{
			this->p1.x = p1.x;
			this->p1.y = p1.y;
			this->p2.x = p2.x;
			this->p2.y = p2.y;

			normal = Utility::Line2Vector(p1, p2).GetAngle() - 90;
			if (normal < 0)
				normal = 360 + normal;
		}

		bool operator==(const Point3D &p1, const Point3D &p2)
		{
			return p1.x == p2.x && p1.y == p2.y && p1.z == p2.z;
		}
		bool operator!=(const Point3D &p1, const Point3D &p2)
		{
			return p1.x == p2.x && p1.y == p2.y && p1.z == p2.z;
		}
		Point2D PointAdd(const Point2D &p1, const Point2D &p2)
		{
			return Point2D(p1.x + p2.x, p1.y + p2.y);
		}
		Point2D operator+(const Point2D &p1, const Point2D &p2)
		{
			return PointAdd(p1, p2);
		}
		Point3D PointAdd(const Point3D &p1, const Point3D &p2)
		{
			return Point3D(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
		}
		Point3D operator+(const Point3D &p1, const Point3D &p2)
		{
			return PointAdd(p1, p2);
		}
		Point2D Transform(Point2D point, Point2D origin, Matrices::MatrixF *transform)
		{
			Vectors::Vector2D vec = Utility::Line2Vector(origin, point);
			vec = vec.Transform(transform);

			return origin + vec.direction;
		}
		Line2D Transform(Line2D line, Point2D origin, Matrices::MatrixF *transform)
		{
			Vectors::Vector2D vec1 = Utility::Line2Vector(origin, line.p1);
			Vectors::Vector2D vec2 = Utility::Line2Vector(origin, line.p2);
			vec1 = vec1.Transform(transform);
			vec2 = vec2.Transform(transform);

			return Line2D(vec1.direction + origin, vec2.direction + origin);
		}
		Point3D RotatePoint(const Point3D &point, Complex::Quaternion &quat)
		{
			Complex::Quaternion retQuat(point);
			retQuat = quat * retQuat * quat.GetInverse();
			return retQuat.GetPoint();
		}

		void PrintProperties(Point2D p)
		{
			printf("X: %f, Y: %f\n", p.x, p.y);
		}
		void PrintProperties(Point3D p)
		{
			printf("X: %.3f, Y: %.3f, Z: %.3f\n", p.x, p.y, p.z);
		}
		void PrintProperties(Line2D l)
		{
			printf("X1: %f, Y1: %f, X2: %f, Y2: %f\n", l.p1.x, l.p1.y, l.p2.x, l.p2.y);
		}
	}
}
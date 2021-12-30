#include "../include/Maths.h"
#include <stdio.h>

namespace MathLib
{
	namespace Primitives
	{
		Line2D::Line2D()
		{
			p1 = Float2();
			p2 = Float2();
			normal = 0;
		}
		Line2D::Line2D(double x1, double y1, double x2, double y2)
		{
			p1.x = x1;
			p1.y = y1;
			p2.x = x2;
			p2.y = y2;
			normal = Utility::Line2Vector(Float2(x1, y1), Float2(x2, y2)).GetAngle() - 90;
		}
		Line2D::Line2D(Float2 p1, Float2 p2)
		{
			this->p1.x = p1.x;
			this->p1.y = p1.y;
			this->p2.x = p2.x;
			this->p2.y = p2.y;

			normal = Utility::Line2Vector(p1, p2).GetAngle() - 90;
			if (normal < 0)
				normal = 360 + normal;
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

		bool operator==(const Int2& f1, const Int2& f2)
		{
			return f1.x == f2.x && f1.y == f2.y;
		}
		bool operator!=(const Int2& f1, const Int2& f2)
		{
			return f1.x != f2.x || f1.y != f2.y;
		}
		bool operator==(const Int3& f1, const Int3& f2)
		{
			return f1.x == f2.x && f1.y == f2.y && f1.z == f2.z;
		}
		bool operator!=(const Int3& f1, const Int3& f2)
		{
			return f1.x != f2.x || f1.y != f2.y || f1.z != f2.z;
		}
		bool operator==(const Int4& f1, const Int4& f2)
		{
			return f1.x == f2.x && f1.y == f2.y && f1.z == f2.z && f1.w == f2.w;
		}
		bool operator!=(const Int4& f1, const Int4& f2)
		{
			return f1.x != f2.x || f1.y != f2.y || f1.z != f2.z || f1.w != f2.w;
		}
		Int2 operator+(const Int2& f1, const Int2& f2)
		{
			return Int2(f1.x + f2.x, f1.y + f2.y);
		}
		Int3 operator+(const Int3& f1, const Int3& f2)
		{
			return Int3(f1.x + f2.x, f1.y + f2.y, f1.z + f2.z);
		}
		Int4 operator+(const Int4& f1, const Int4& f2)
		{
			return Int4(f1.x + f2.x, f1.y + f2.y, f1.z + f2.z, f1.w + f2.w);
		}
		Int2 operator-(const Int2& f1, const Int2& f2)
		{
			return Int2(f1.x - f2.x, f1.y - f2.y);
		}
		Int3 operator-(const Int3& f1, const Int3& f2)
		{
			return Int3(f1.x - f2.x, f1.y - f2.y, f1.z - f2.z);
		}
		Int4 operator-(const Int4& f1, const Int4& f2)
		{
			return Int4(f1.x - f2.x, f1.y - f2.y, f1.z - f2.z, f1.w - f2.w);
		}

		Float2 Transform(Float2 point, Float2 origin, const Matrices::MatrixF& transform)
		{
			Vectors::Vector2D vec = Utility::Line2Vector(origin, point);
			vec = vec.Transform(transform);

			return origin + vec.direction;
		}
		Line2D Transform(Line2D line, Float2 origin, const Matrices::MatrixF& transform)
		{
			Vectors::Vector2D vec1 = Utility::Line2Vector(origin, line.p1);
			Vectors::Vector2D vec2 = Utility::Line2Vector(origin, line.p2);
			vec1 = vec1.Transform(transform);
			vec2 = vec2.Transform(transform);

			return Line2D(vec1.direction + origin, vec2.direction + origin);
		}
		Float3 RotatePoint(const Float3& point, Complex::Quaternion& quat)
		{
			Complex::Quaternion retQuat(point);
			retQuat = quat * retQuat * quat.GetInverse();
			return retQuat.GetPoint();
		}

		void PrintProperties(Float2 p)
		{
			printf("X: %f, Y: %f\n", p.x, p.y);
		}
		void PrintProperties(Float3 p)
		{
			printf("X: %.3f, Y: %.3f, Z: %.3f\n", p.x, p.y, p.z);
		}
		void PrintProperties(Line2D l)
		{
			printf("X1: %f, Y1: %f, X2: %f, Y2: %f\n", l.p1.x, l.p1.y, l.p2.x, l.p2.y);
		}
	}
}
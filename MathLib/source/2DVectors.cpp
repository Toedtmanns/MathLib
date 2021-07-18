#include "../include/Maths.h"

namespace MathLib
{
	namespace Vectors
	{
		Vector2D::Vector2D() : direction({0, 0})
		{

		}
		Vector2D::Vector2D(Primitives::Point2D dir) : direction(dir)
		{

		}
		Vector2D::Vector2D(double angle, double length) : direction({0, 0})
		{
			direction.x += sin(angle * PI / 180) * length;
			direction.y += cos(angle * PI / 180) * length;
		}
		Vector2D Vector2D::Transform(Matrices::MatrixF **transformation)
		{
			if (!Matrices::MatrixIsSquare(transformation, 2))
				throw std::invalid_argument("Matrix is not 2x2!");

			double dirX, dirY;

			dirX = direction.x * (*transformation)->GetRow(0).GetAt(0) + direction.y * (*transformation)->GetRow(1).GetAt(0);
			dirY = direction.x * (*transformation)->GetRow(0).GetAt(1) + direction.y * (*transformation)->GetRow(1).GetAt(1);

			return Vector2D(Primitives::Point2D(dirX, dirY));
		}
		Vector2D Vector2D::Scale(double s)
		{
			return Vector2D(Primitives::Point2D(direction.x * s, direction.y * s));
		}
		Vector2D Vector2D::Scale(double sX, double sY)
		{
			return Vector2D(Primitives::Point2D(direction.x * sX, direction.y * sY));
		}
		Vector2D Vector2D::Rotate(double angle)
		{
			if (fmod(angle, 180) == 0)
				return Vector2D(Primitives::Point2D(-direction.x, -direction.y));

			angle = -angle;
			Primitives::Point2D p1 = Primitives::Point2D(
				sinf(Utility::Deg2Rad(angle + 90.0)),
				cosf(Utility::Deg2Rad(angle + 90.0))
			);
			Primitives::Point2D p2 = Primitives::Point2D(
				sinf(Utility::Deg2Rad(angle)),
				cosf(Utility::Deg2Rad(angle))
			);
			Matrices::MatrixF *transform = Matrices::TransformF2x2(p1, p2);
			return Transform(&transform);
		}
		double Vector2D::GetAngle()
		{
			bool neg = false;
			if (direction.x < 0)
			{
				neg = true;
				direction.x *= -1;
			}

			double hyp = this->GetLen();
			double angle = Utility::Rad2Deg(acosf(direction.y / hyp));

			if (neg)
				angle = 360 - angle;

			return angle;
		}
		double Vector2D::GetLen()
		{
			return sqrtf(powf(direction.x, 2) + powf(direction.y, 2));
		}
		Primitives::Point2D Vector2D::ApplyVector(Primitives::Point2D p)
		{
			p.x += direction.x;
			p.y += direction.y;

			return p;
		}

		Vector2D VectorAdd(Vector2D v1, Vector2D v2)
		{
			v1.direction.x += v2.direction.x;
			v1.direction.y += v2.direction.y;
			return v1;
		}
		Vector2D VectorSub(Vector2D v1, Vector2D v2)
		{
			v2 = v2.Rotate(180);
			return VectorAdd(v1, v2);
		}
		double VectorGetAngleDifference(Vector2D v1, Vector2D v2)
		{
			double a1 = v1.GetAngle();
			double a2 = v2.GetAngle();

			double angle = (double)fmaxf(a1, a2) - (double)fminf(a1, a2);
			if (angle > 180)
				angle = 360 - angle;

			return angle;
		}
		double VectorDotProduct(Vector2D v1, Vector2D v2)
		{
			return v1.direction.x * v2.direction.x + v1.direction.y * v2.direction.y;
		}
		Vector3D VectorCrossProduct(Vector2D v1, Vector2D v2)
		{
			Matrices::MatrixF *vectors = Matrices::TransformF2x2(v1.direction, v2.direction);
			double zLength = Matrices::MatrixGetDet(&vectors);
			Vector3D cProduct = Vector3D(Primitives::Point3D(0, 0, zLength));

			return cProduct;
		}
		double VectorGetDeterminant(Vector2D v1, Vector2D v2)
		{
			Matrices::MatrixF *vectors = Matrices::TransformF2x2(v1.direction, v2.direction);
			return Matrices::MatrixGetDet(&vectors);
		}
		void PrintProperties(Vector2D v)
		{
			printf("X: %.3f, Y: %.3f\n", v.direction.x, v.direction.y);
		}
	}
}
#include "../include/Maths.h"

namespace MathLib
{
	namespace Vectors
	{
		Vector3D::Vector3D()
			: direction(Primitives::Point3D(0, 0, 0))
		{

		}
		Vector3D::Vector3D(Primitives::Point3D dir)
			: direction(dir)
		{

		}
		Vector3D Vector3D::Transform(Matrices::MatrixF *mat)
		{
			if (!Matrices::MatrixIsSquare(mat, 3))
				throw std::invalid_argument("Matrix is not 3x3!");

			double dirX, dirY, dirZ;

			dirX = direction.x * mat->GetRow(0).GetAt(0) + direction.y * mat->GetRow(1).GetAt(0) + direction.z * mat->GetRow(2).GetAt(0);
			dirY = direction.x * mat->GetRow(0).GetAt(1) + direction.y * mat->GetRow(1).GetAt(1) + direction.z * mat->GetRow(2).GetAt(1);
			dirZ = direction.x * mat->GetRow(0).GetAt(2) + direction.y * mat->GetRow(1).GetAt(2) + direction.z * mat->GetRow(2).GetAt(2);

			return Vector3D(Primitives::Point3D(dirX, dirY, dirZ));
		}
		Vector3D Vector3D::Scale(double s)
		{
			return Vector3D(Primitives::Point3D(direction.x * s, direction.y * s, direction.z * s));
		}
		Vector3D Vector3D::Scale(double sX, double sY, double sZ)
		{
			return Vector3D(Primitives::Point3D(direction.x * sX, direction.y * sY, direction.z * sZ));
		}
		Vector3D Vector3D::Rotate(double angle, int axis)
		{
			Vector3D retVec = *this;
			Vector2D calcVec = Vector2D();

			switch (axis)
			{
			case 0:
				calcVec.direction = Primitives::Point2D(direction.y, direction.z);
				calcVec = calcVec.Rotate(angle);
				retVec.direction = Primitives::Point3D(direction.x, calcVec.direction.x, calcVec.direction.y);
				break;
			case 1:
				calcVec.direction = Primitives::Point2D(direction.x, direction.z);
				calcVec = calcVec.Rotate(angle);
				retVec.direction = Primitives::Point3D(calcVec.direction.x, direction.y, calcVec.direction.y);
				break;
			case 2:
				calcVec.direction = Primitives::Point2D(direction.x, direction.y);
				calcVec = calcVec.Rotate(angle);
				retVec.direction = Primitives::Point3D(calcVec.direction.x, calcVec.direction.y, direction.z);
				break;
			default:
				throw std::invalid_argument("There are only 3 axes!");
			}

			return retVec;
		}
		/*Vector3D Vector3D::Rotate(double angle, Vector3D *normal)
		{

		}*/
		double Vector3D::GetAngle(int axis)
		{
			Vector2D calcVec = Vector2D();
			double angle;

			switch (axis)
			{
			case 0:
				calcVec.direction = Primitives::Point2D(direction.y, direction.z);
				angle = calcVec.GetAngle();
				break;
			case 1:
				calcVec.direction = Primitives::Point2D(direction.x, direction.z);
				angle = calcVec.GetAngle();
				break;
			case 2:
				calcVec.direction = Primitives::Point2D(direction.x, direction.y);
				angle = calcVec.GetAngle();
				break;
			default:
				throw std::invalid_argument("There are only 3 axes!");
			}

			return angle;
		}
		/*double Vector3D::GetAngle(Vector3D *normal)
		{

		}*/
		double Vector3D::GetLen()
		{
			return sqrt(powf(sqrtf(powf(direction.x, 2) + powf(direction.y, 2)), 2) + powf(direction.z, 2));
		}
		void PrintProperties(Vector3D v)
		{
			printf("X: %.3f, Y: %.3f, Z: %.3f\n", v.direction.x, v.direction.y, v.direction.z);
		}
	}
}
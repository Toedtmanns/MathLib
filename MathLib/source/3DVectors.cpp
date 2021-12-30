#include "../include/Maths.h"
#include <stdexcept>
#include <cmath>

namespace MathLib
{
	namespace Vectors
	{
		Vector3D::Vector3D()
			: direction(Primitives::Float3(0, 0, 0))
		{

		}
		Vector3D::Vector3D(Primitives::Float3 dir)
			: direction(dir)
		{

		}
		Vector3D Vector3D::Transform(const Matrices::MatrixF& mat)
		{
			if (!Matrices::MatrixIsSquare(mat, 3))
				throw std::invalid_argument("Matrix is not 3x3!");

			double dirX, dirY, dirZ;

			dirX = direction.x * mat.GetNum(0, 0) + direction.y * mat.GetNum(0, 1) + direction.z * mat.GetNum(0, 2);
			dirY = direction.x * mat.GetNum(1, 0) + direction.y * mat.GetNum(1, 1) + direction.z * mat.GetNum(1, 2);
			dirZ = direction.x * mat.GetNum(1, 0) + direction.y * mat.GetNum(2, 1) + direction.z * mat.GetNum(2, 2);

			return Vector3D(Primitives::Float3(dirX, dirY, dirZ));
		}
		Vector3D Vector3D::Scale(double s)
		{
			return Vector3D(Primitives::Float3(direction.x * s, direction.y * s, direction.z * s));
		}
		Vector3D Vector3D::Scale(double sX, double sY, double sZ)
		{
			return Vector3D(Primitives::Float3(direction.x * sX, direction.y * sY, direction.z * sZ));
		}
		Vector3D Vector3D::Rotate(double angle, int axis)
		{
			Vector3D retVec = *this;
			Vector2D calcVec = Vector2D();

			switch (axis)
			{
			case 0:
				calcVec.direction = Primitives::Float2(direction.y, direction.z);
				calcVec = calcVec.Rotate(angle);
				retVec.direction = Primitives::Float3(direction.x, calcVec.direction.x, calcVec.direction.y);
				break;
			case 1:
				calcVec.direction = Primitives::Float2(direction.x, direction.z);
				calcVec = calcVec.Rotate(angle);
				retVec.direction = Primitives::Float3(calcVec.direction.x, direction.y, calcVec.direction.y);
				break;
			case 2:
				calcVec.direction = Primitives::Float2(direction.x, direction.y);
				calcVec = calcVec.Rotate(angle);
				retVec.direction = Primitives::Float3(calcVec.direction.x, calcVec.direction.y, direction.z);
				break;
			default:
				throw std::invalid_argument("There are only 3 axes!");
			}

			return retVec;
		}
		/*Vector3D Vector3D::Rotate(double angle, Vector3D *normal)
		{

		}*/
		double Vector3D::GetAngle(int axis) const
		{
			Vector2D calcVec = Vector2D();
			double angle;

			switch (axis)
			{
			case 0:
				calcVec.direction = Primitives::Float2(direction.y, direction.z);
				angle = calcVec.GetAngle();
				break;
			case 1:
				calcVec.direction = Primitives::Float2(direction.x, direction.z);
				angle = calcVec.GetAngle();
				break;
			case 2:
				calcVec.direction = Primitives::Float2(direction.x, direction.y);
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
		double Vector3D::GetLen() const
		{
			return sqrt(powf(sqrtf(powf(direction.x, 2) + powf(direction.y, 2)), 2) + powf(direction.z, 2));
		}
		Vector3D Vector3D::SetLen(double len)
		{
			double vecLen(GetLen());
			Vector3D retVec(*this);
			return retVec.Scale(vecLen * (len / vecLen));
		}
		void PrintProperties(Vector3D v)
		{
			printf("X: %.3f, Y: %.3f, Z: %.3f\n", v.direction.x, v.direction.y, v.direction.z);
		}
		Vector3D GetRelativeVec(const Vector3D& vec1, const Vector3D& vec2)
		{
			Vector3D retVec = Vector3D();
			retVec.direction.x = vec2.direction.x - vec1.direction.x;
			retVec.direction.y = vec2.direction.y - vec1.direction.y;
			retVec.direction.z = vec2.direction.z - (vec1.direction.z - 1);
			return retVec;
		}
	}
}
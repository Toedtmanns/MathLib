#include "../include/Maths.h"
#include <stdexcept>
#include <cmath>

namespace MathLib
{
	namespace Vectors
	{
		// 2D Vectors

		Vector2D::Vector2D() : direction({0, 0})
		{

		}
		Vector2D::Vector2D(const Primitives::Float2& dir)
			: direction(dir)
		{

		}
		Vector2D::Vector2D(const double& angle, const double& length)
			: direction({0, 0})
		{
			direction.x += cos(angle * PI / 180) * length;
			direction.y += sin(angle * PI / 180) * length;
		}
		Vector2D::Vector2D(const Primitives::Line2D& line)
			: direction({line.p2.x - line.p1.x, line.p2.y - line.p1.y})
		{
			
		}
		void Vector2D::Transform(const Matrices::MatrixF& transformMat)
		{
			if (!Matrices::MatrixIsSquare(transformMat, 2))
				throw std::invalid_argument("Matrix is not 2x2!");

			direction = direction * transformMat;
		}
		void Vector2D::Scale(const double& scale)
		{
			direction.x *= scale;
			direction.y *= scale;
		}
		void Vector2D::Scale(const double& scaleX, const double& scaleY)
		{
			direction.x *= scaleX;
			direction.y *= scaleY;
		}
		void Vector2D::SetScale(const double& scale)
		{
			double factor = sqrt(pow(direction.x, 2) + pow(direction.y, 2));
			direction = direction * (1.0 / factor);
		}
		void Vector2D::Rotate(double angle)
		{
			Matrices::MatrixF rotMat = Matrices::MatrixF(2);
			rotMat[0][0] = sin(Utility::Deg2Rad(angle + 90));
			rotMat[0][1] = cos(Utility::Deg2Rad(angle + 90));
			rotMat[1][0] = sin(Utility::Deg2Rad(angle));
			rotMat[1][1] = cos(Utility::Deg2Rad(angle));

			direction = direction * rotMat;
		}
		Vector2D Vector2D::operator+(const Vector2D& other) const
		{
			return Vector2D(Primitives::Float2(direction.x + other.direction.x, direction.y + other.direction.y));
		}
		Vector2D Vector2D::operator-(const Vector2D& other) const
		{
			return Vector2D({direction.x - other.direction.x, direction.y - other.direction.y});
		}
		double Vector2D::operator*(const Vector2D& other) const
		{
			return direction.x * other.direction.x + direction.y * other.direction.y;
		}
		double Vector2D::GetAngle() const
		{
			double hyp = this->GetLen();
			double angle = Utility::Rad2Deg(acos(direction.x / hyp));

			if (direction.y > 0)
				angle = -angle;

			return angle;
		}
		double Vector2D::GetLen() const
		{
			return sqrt(pow(direction.x, 2) + pow(direction.y, 2));
		}
		Matrices::MatrixF Vector2D::GetRowVector() const
		{
			Matrices::MatrixF retMat(2, 1);
			retMat.SetNum(0, 0, direction.x);
			retMat.SetNum(1, 0, direction.y);

			return retMat;
		}
		Matrices::MatrixF Vector2D::GetColVector() const
		{
			Matrices::MatrixF retMat(1, 2);
			retMat.SetNum(0, 0, direction.x);
			retMat.SetNum(0, 1, direction.y);

			return retMat;
		}
		Primitives::Float2 Vector2D::TransformPoint(Primitives::Float2 point) const
		{
			point.x += direction.x;
			point.y += direction.y;

			return point;
		}

		// 3D Vectors

		Vector3D::Vector3D()
			: direction({0, 0, 0})
		{

		}
		Vector3D::Vector3D(const Primitives::Float3& dir)
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

			direction.x = dirX;
			direction.y = dirY;
			direction.z = dirZ;
		}
		Vector3D Vector3D::Scale(const double& s)
		{
			return Vector3D(Primitives::Float3(direction.x * s, direction.y * s, direction.z * s));
		}
		Vector3D Vector3D::Scale(const double& sX, const double& sY, const double& sZ)
		{
			return Vector3D(Primitives::Float3(direction.x * sX, direction.y * sY, direction.z * sZ));
		}
		void Vector3D::Rotate(const double& angle, const unsigned int& axis)
		{
			Vector2D calcVec = Vector2D();

			switch (axis % 3)
			{
			case 0:
				calcVec.direction = Primitives::Float2(direction.y, direction.z);
				calcVec.Rotate(angle);
				direction = Primitives::Float3(direction.x, calcVec.direction.x, calcVec.direction.y);
				break;
			case 1:
				calcVec.direction = Primitives::Float2(direction.x, direction.z);
				calcVec.Rotate(angle);
				direction = Primitives::Float3(calcVec.direction.x, direction.y, calcVec.direction.y);
				break;
			case 2:
				calcVec.direction = Primitives::Float2(direction.x, direction.y);
				calcVec.Rotate(angle);
				direction = Primitives::Float3(calcVec.direction.x, calcVec.direction.y, direction.z);
				break;
			}
		}
		/*Vector3D Vector3D::Rotate(double angle, Vector3D *normal)
		{

		}*/
		Vector3D Vector3D::SetLen(const double& len)
		{
			double vecLen(GetLen());
			Vector3D retVec(*this);

			return retVec.Scale(vecLen * (len / vecLen));
		}
		Vector3D Vector3D::operator+(const Vector3D& other) const
		{
			return Vector3D({
				direction.x + other.direction.x,
				direction.y + other.direction.y,
				direction.z + other.direction.z
				});
		}
		Vector3D Vector3D::operator-(const Vector3D& other) const
		{
			return Vector3D({
				direction.x - other.direction.x,
				direction.y - other.direction.y,
				direction.z - other.direction.z
				});
		}
		double Vector3D::operator*(const Vector3D& other) const
		{
			return direction.x * other.direction.x + direction.y * other.direction.y + direction.z * other.direction.z;
		}
		double Vector3D::GetAngle(const int& axis) const
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
			return sqrt(pow(sqrt(pow(direction.x, 2) + pow(direction.y, 2)), 2) + pow(direction.z, 2));
		}
		Matrices::MatrixF Vector3D::GetRowVector() const
		{
			Matrices::MatrixF retMat(3, 1);
			retMat.SetNum(0, 0, direction.x);
			retMat.SetNum(1, 0, direction.y);
			retMat.SetNum(2, 0, direction.z);

			return retMat;
		}
		Matrices::MatrixF Vector3D::GetColVector() const
		{
			Matrices::MatrixF retMat(1, 3);
			retMat.SetNum(0, 0, direction.x);
			retMat.SetNum(0, 1, direction.y);
			retMat.SetNum(0, 2, direction.z);

			return retMat;
		}
		Primitives::Float3 Vector3D::TransformPoint(Primitives::Float3 point)
		{
			point.x += direction.x;
			point.y += direction.y;
			point.z += direction.z;

			return point;
		}

		Vector3D GetRelativeVec(const Vector3D& vec1, const Vector3D& vec2)
		{
			Vector3D retVec = Vector3D();
			retVec.direction.x = vec2.direction.x - vec1.direction.x;
			retVec.direction.y = vec2.direction.y - vec1.direction.y;
			retVec.direction.z = vec2.direction.z - (vec1.direction.z - 1);

			return retVec;
		}
		double VectorGetAngleDifference(const Vector2D& v1, const Vector2D& v2)
		{
			double a1 = v1.GetAngle();
			double a2 = v2.GetAngle();

			double angle = (double) fmax(a1, a2) - (double) fmin(a1, a2);
			if (angle > 180)
				angle = 360 - angle;

			return angle;
		}
		double VectorDotProduct(const Vector2D& v1, const Vector2D& v2)
		{
			return v1.direction.x * v2.direction.x + v1.direction.y * v2.direction.y;
		}
		Vector3D VectorCrossProduct(const Vector2D& v1, const Vector2D& v2)
		{
			Matrices::MatrixF vectors = Matrices::TransformF2x2(v1.direction, v2.direction);
			double zLength = Matrices::MatrixGetDet(vectors);
			Vector3D cProduct = Vector3D(Primitives::Float3(0, 0, zLength));

			return cProduct;
		}
		double VectorGetDeterminant(const Vector2D& v1, const Vector2D& v2)
		{
			Matrices::MatrixF vectors = Matrices::TransformF2x2(v1.direction, v2.direction);
			return Matrices::MatrixGetDet(vectors);
		}
		void PrintProperties(const Vector2D& vector)
		{
			printf("X: %.3f, Y: %.3f\n", vector.direction.x, vector.direction.y);
		}
		void PrintProperties(const Vector3D& vector)
		{
			printf("X: %.3f, Y: %.3f, Z: %.3f\n", vector.direction.x, vector.direction.y, vector.direction.z);
		}
	}
}
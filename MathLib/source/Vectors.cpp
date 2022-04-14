#include "../include/Maths.h"
#include <stdexcept>
#include <cmath>

namespace MathLib
{
	// 2D Vectors

	Vector2D::Vector2D() : direction({0, 0})
	{

	}
	Vector2D::Vector2D(const Float2& dir)
		: direction(dir)
	{

	}
	Vector2D::Vector2D(const double& x, const double& y)
		: direction(x, y)
	{

	}
	Vector2D::Vector2D(const Line2D& line)
		: direction({line.p2.x - line.p1.x, line.p2.y - line.p1.y})
	{
		
	}
	Vector2D::Vector2D(const Float2& p1, const Float2& p2)
		: direction({p2.x - p1.x, p2.y - p1.y})
	{

	}
	void Vector2D::Transform(const MatrixF& transformMat)
	{
		if (!MatrixIsSquare(transformMat, 2))
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
		MatrixF rotMat = MatrixF(2);
		rotMat[0][0] = sin(Deg2Rad(angle + 90));
		rotMat[0][1] = cos(Deg2Rad(angle + 90));
		rotMat[1][0] = sin(Deg2Rad(angle));
		rotMat[1][1] = cos(Deg2Rad(angle));

		direction = direction * rotMat;
	}
	void Vector2D::Rot90R()
	{
		double tmp = direction.x;
		direction.x = direction.y;
		direction.y = -tmp;
	}
	void Vector2D::Rot90L()
	{
		double tmp = direction.x;
		direction.x = -direction.y;
		direction.y = tmp;
	}
	Vector2D Vector2D::TripleProduct(const Vector2D& other) const
	{
		return VectorTripleProduct(*this, other);
	}
	Vector2D Vector2D::operator+(const Vector2D& other) const
	{
		return Vector2D(Float2(direction.x + other.direction.x, direction.y + other.direction.y));
	}
	Vector2D Vector2D::operator-(const Vector2D& other) const
	{
		return Vector2D({direction.x - other.direction.x, direction.y - other.direction.y});
	}
	double Vector2D::operator*(const Vector2D& other) const
	{
		return direction.x * other.direction.x + direction.y * other.direction.y;
	}
	Vector2D Vector2D::operator*(const double& number) const
	{
		return {
			direction.x * number,
			direction.y * number
		};
	}
	void Vector2D::operator*=(const double& number)
	{
		direction.x *= number;
		direction.y *= number;
	}
	Vector2D Vector2D::operator-() const
	{
		return {-direction.x, -direction.y};
	}
	double Vector2D::GetAngle() const
	{
		double hyp = this->GetLen();
		double angle = fmod(180.001 + Rad2Deg(acos(direction.y / hyp)), 180.001);

		if (direction.x < 0)
			angle = -angle;

		return angle;
	}
	double Vector2D::GetLen() const
	{
		return sqrt(pow(direction.x, 2) + pow(direction.y, 2));
	}
	MatrixF Vector2D::GetRowVector() const
	{
		MatrixF retMat(2, 1);
		retMat.SetNum(0, 0, direction.x);
		retMat.SetNum(1, 0, direction.y);

		return retMat;
	}
	MatrixF Vector2D::GetColVector() const
	{
		MatrixF retMat(1, 2);
		retMat.SetNum(0, 0, direction.x);
		retMat.SetNum(0, 1, direction.y);

		return retMat;
	}
	Float2 Vector2D::TransformPoint(Float2 point) const
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
	Vector3D::Vector3D(const Float3& dir)
		: direction(dir)
	{

	}
	Vector3D::Vector3D(const double& x, const double& y, const double& z)
		: direction(x, y, z)
	{

	}
	Vector3D::Vector3D(const Vector2D& vec, const double& z)
		: direction(vec.direction.x, vec.direction.y, z)
	{

	}
	void Vector3D::Transform(const MatrixF& mat)
	{
		if (!MatrixIsSquare(mat, 3))
			throw std::invalid_argument("Matrix is not 3x3!");

		double dirX, dirY, dirZ;

		dirX = direction.x * mat.GetNum(0, 0) + direction.y * mat.GetNum(0, 1) + direction.z * mat.GetNum(0, 2);
		dirY = direction.x * mat.GetNum(1, 0) + direction.y * mat.GetNum(1, 1) + direction.z * mat.GetNum(1, 2);
		dirZ = direction.x * mat.GetNum(1, 0) + direction.y * mat.GetNum(2, 1) + direction.z * mat.GetNum(2, 2);

		direction.x = dirX;
		direction.y = dirY;
		direction.z = dirZ;
	}
	void Vector3D::Scale(const double& s)
	{
		direction.x *= s;
		direction.y *= s;
		direction.z *= s;
	}
	void Vector3D::Scale(const double& sX, const double& sY, const double& sZ)
	{
		direction.x *= sX;
		direction.y *= sY;
		direction.z *= sZ;
	}
	void Vector3D::Rotate(const double& angle, const unsigned int& axis)
	{
		Vector2D calcVec = Vector2D();

		switch (axis % 3)
		{
		case 0:
			calcVec.direction = Float2(direction.y, direction.z);
			calcVec.Rotate(angle);
			direction = Float3(direction.x, calcVec.direction.x, calcVec.direction.y);
			break;
		case 1:
			calcVec.direction = Float2(direction.x, direction.z);
			calcVec.Rotate(angle);
			direction = Float3(calcVec.direction.x, direction.y, calcVec.direction.y);
			break;
		case 2:
			calcVec.direction = Float2(direction.x, direction.y);
			calcVec.Rotate(angle);
			direction = Float3(calcVec.direction.x, calcVec.direction.y, direction.z);
			break;
		}
	}
	/*Vector3D Vector3D::Rotate(double angle, Vector3D *normal)
	{

	}*/
	void Vector3D::SetLen(const double& len)
	{
		double vecLen = GetLen();

		Scale(vecLen * (len / vecLen));
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
	Vector3D Vector3D::operator*(const double& number) const
	{
		return Vector3D(direction.x * number, direction.y * number, direction.z * number);
	}
	void Vector3D::operator*=(const double& number)
	{
		direction.x *= number;
		direction.y *= number;
		direction.z *= number;
	}
	Vector3D Vector3D::operator-() const
	{
		return Vector3D(-direction.x, -direction.y, -direction.z);
	}
	double Vector3D::GetAngle(const int& axis) const
	{
		Vector2D calcVec = Vector2D();
		double angle;

		switch (axis)
		{
		case 0:
			calcVec.direction = Float2(direction.y, direction.z);
			angle = calcVec.GetAngle();
			break;
		case 1:
			calcVec.direction = Float2(direction.x, direction.z);
			angle = calcVec.GetAngle();
			break;
		case 2:
			calcVec.direction = Float2(direction.x, direction.y);
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
	MatrixF Vector3D::GetRowVector() const
	{
		MatrixF retMat(3, 1);
		retMat.SetNum(0, 0, direction.x);
		retMat.SetNum(1, 0, direction.y);
		retMat.SetNum(2, 0, direction.z);

		return retMat;
	}
	MatrixF Vector3D::GetColVector() const
	{
		MatrixF retMat(1, 3);
		retMat.SetNum(0, 0, direction.x);
		retMat.SetNum(0, 1, direction.y);
		retMat.SetNum(0, 2, direction.z);

		return retMat;
	}
	Float3 Vector3D::TransformPoint(Float3 point)
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
		MatrixF vectors = TransformF2x2(v1.direction, v2.direction);
		double zLength = MatrixGetDet(vectors);
		Vector3D cProduct = Vector3D(Float3(0, 0, zLength));

		return cProduct;
	}
	Vector3D VectorCrossProduct(const Vector3D& v1, const Vector3D& v2)
	{
		return Vector3D{
			v1.direction.y * v2.direction.z - v1.direction.z * v2.direction.y,
			v1.direction.z * v2.direction.x - v1.direction.x * v2.direction.z,
			v1.direction.x * v2.direction.y - v1.direction.y * v2.direction.x
		};
	}
	Vector2D VectorTripleProduct(const Vector2D& v1, const Vector2D& v2)
	{
		Vector3D tpVec = VectorCrossProduct(VectorCrossProduct(Vector3D(v1), Vector3D(v2)), Vector3D(v1));
		return Vector2D(tpVec.direction.x, tpVec.direction.y);
	}
	Vector3D VectorTripleProduct(const Vector3D& v1, const Vector3D& v2)
	{
		return VectorCrossProduct(VectorCrossProduct(Vector3D(v1), Vector3D(v2)), Vector3D(v1));
	}
	double VectorGetDeterminant(const Vector2D& v1, const Vector2D& v2)
	{
		MatrixF vectors = TransformF2x2(v1.direction, v2.direction);
		return MatrixGetDet(vectors);
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
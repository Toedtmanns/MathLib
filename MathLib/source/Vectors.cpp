#include "../include/Maths.hpp"
#include <stdexcept>
#include <cmath>

namespace MathLib
{
	// 2D Vectors

	Vector2D::Vector2D() 
		: Float2(0, 0)
	{

	}
	Vector2D::Vector2D(const Float2& dir)
		: Float2(dir)
	{

	}
	Vector2D::Vector2D(const double& x, const double& y)
		: Float2(x, y)
	{

	}
	Vector2D::Vector2D(const Line2D& line)
		: Float2({line.p2.x - line.p1.x, line.p2.y - line.p1.y})
	{
		
	}
	Vector2D::Vector2D(const Float2& p1, const Float2& p2)
		: Float2({p2.x - p1.x, p2.y - p1.y})
	{

	}
	const Vector2D& Vector2D::Transform(const MatrixF& transformMat)
	{
		if (!MatrixIsSquare(transformMat, 2))
			throw std::invalid_argument("Matrix is not 2x2!");

		*this = *this * transformMat;
		return *this;
	}
	const Vector2D& Vector2D::Scale(const double& scale)
	{
		x *= scale;
		y *= scale;
		return *this;
	}
	const Vector2D& Vector2D::Scale(const double& scaleX, const double& scaleY)
	{
		x *= scaleX;
		y *= scaleY;
		return *this;
	}
	const Vector2D& Vector2D::SetScale(const double& scale)
	{
		double factor = sqrt(pow(x, 2) + pow(y, 2));
		*this = *this * (1.0 / factor);
		return *this;
	}
	const Vector2D& Vector2D::Rotate(double angle)
	{
		MatrixF rotMat = MatrixF(2);
		rotMat[0][0] = sin(Deg2Rad(angle + 90));
		rotMat[0][1] = cos(Deg2Rad(angle + 90));
		rotMat[1][0] = sin(Deg2Rad(angle));
		rotMat[1][1] = cos(Deg2Rad(angle));

		*this = *this * rotMat;
		return *this;
	}
	const Vector2D& Vector2D::Rot90R()
	{
		double tmp = x;
		x = y;
		y = -tmp;
		return *this;
	}
	const Vector2D& Vector2D::Rot90L()
	{
		double tmp = x;
		x = -y;
		y = tmp;
		return *this;
	}
	Vector2D Vector2D::TripleProduct(const Vector2D& other) const
	{
		return VectorTripleProduct(*this, other);
	}
	Vector2D Vector2D::operator+(const Vector2D& other) const
	{
		return Vector2D(Float2(x + other.x, y + other.y));
	}
	Vector2D Vector2D::operator-(const Vector2D& other) const
	{
		return Vector2D({x - other.x, y - other.y});
	}
	double Vector2D::operator*(const Vector2D& other) const
	{
		return x * other.x + y * other.y;
	}
	Vector2D Vector2D::operator*(const double& number) const
	{
		return {
			x * number,
			y * number
		};
	}
	Vector2D Vector2D::operator*(const MatrixF& matrix) const
	{
		MatrixF result = GetColVector() * matrix;
		return {result[0][0], result[0][1]};
	}
	void Vector2D::operator*=(const double& number)
	{
		x *= number;
		y *= number;
	}
	Vector2D Vector2D::operator-() const
	{
		return {-x, -y};
	}
	double Vector2D::GetAngle() const
	{
		double hyp = this->GetLen();
		double angle = fmod(180.001 + Rad2Deg(acos(y / hyp)), 180.001);

		if (x < 0)
			angle = -angle;

		return angle;
	}
	double Vector2D::GetLen() const
	{
		return sqrt(pow(x, 2) + pow(y, 2));
	}
	MatrixF Vector2D::GetRowVector() const
	{
		MatrixF retMat(2, 1);
		retMat.SetNum(0, 0, x);
		retMat.SetNum(1, 0, y);

		return retMat;
	}
	MatrixF Vector2D::GetColVector() const
	{
		MatrixF retMat(1, 2);
		retMat.SetNum(0, 0, x);
		retMat.SetNum(0, 1, y);

		return retMat;
	}
	Float2 Vector2D::TransformPoint(Float2 point) const
	{
		point.x += x;
		point.y += y;

		return point;
	}

	// 3D Vectors

	Vector3D::Vector3D()
		: Float3(0, 0, 0)
	{

	}
	Vector3D::Vector3D(const Float3& dir)
		: Float3(dir)
	{

	}
	Vector3D::Vector3D(const double& x, const double& y, const double& z)
		: Float3(x, y, z)
	{

	}
	Vector3D::Vector3D(const Vector2D& vec, const double& z)
		: Float3(vec.x, vec.y, z)
	{

	}
	Vector3D::Vector3D(const Float3& p1, const Float3& p2)
		: Float3(p2 - p1)
	{

	}
	const Vector3D& Vector3D::Transform(const MatrixF& mat)
	{
		MatrixF result = GetColVector() * mat;
		x = result[0][0];
		y = result[0][1];
		z = result[0][2];
		return *this;
	}
	const Vector3D& Vector3D::Scale(const double& s)
	{
		x *= s;
		y *= s;
		z *= s;
		return *this;
	}
	const Vector3D& Vector3D::Scale(const double& sX, const double& sY, const double& sZ)
	{
		x *= sX;
		y *= sY;
		z *= sZ;
		return *this;
	}
	const Vector3D& Vector3D::SetScale(const double& len)
	{
		double vecLen = GetLen();

		Scale(vecLen * (len / vecLen));
		return *this;
	}
	const Vector3D& Vector3D::Rotate(const double& angle, const unsigned int& axis)
	{
		Vector2D calcVec = Vector2D();

		switch (axis % 3)
		{
		case 0:
			calcVec = Vector2D(y, z);
			calcVec.Rotate(angle);
			*this = Vector3D(x, calcVec.x, calcVec.y);
			break;
		case 1:
			calcVec = Vector2D(x, z);
			calcVec.Rotate(angle);
			*this = Vector3D(calcVec.x, y, calcVec.y);
			break;
		case 2:
			calcVec = Vector2D(x, y);
			calcVec.Rotate(angle);
			*this = Vector3D(calcVec.x, calcVec.y, z);
			break;
		}
		return *this;
	}
	/*Vector3D Vector3D::Rotate(double angle, Vector3D *normal)
	{

	}*/
	Vector3D Vector3D::operator+(const Vector3D& other) const
	{
		return Vector3D({
			x + other.x,
			y + other.y,
			z + other.z
			});
	}
	Vector3D Vector3D::operator-(const Vector3D& other) const
	{
		return Vector3D({
			x - other.x,
			y - other.y,
			z - other.z
			});
	}
	double Vector3D::operator*(const Vector3D& other) const
	{
		return x * other.x + y * other.y + z * other.z;
	}
	Vector3D Vector3D::operator*(const double& number) const
	{
		return Vector3D(x * number, y * number, z * number);
	}
	void Vector3D::operator*=(const double& number)
	{
		x *= number;
		y *= number;
		z *= number;
	}
	Vector3D Vector3D::operator-() const
	{
		return Vector3D(-x, -y, -z);
	}
	double Vector3D::GetAngle(const int& axis) const
	{
		Vector2D calcVec = Vector2D();
		double angle;

		switch (axis)
		{
		case 0:
			calcVec = Vector2D(y, z);
			angle = calcVec.GetAngle();
			break;
		case 1:
			calcVec = Vector2D(x, z);
			angle = calcVec.GetAngle();
			break;
		case 2:
			calcVec = Vector2D(x, y);
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
		return sqrt(pow(sqrt(pow(x, 2) + pow(y, 2)), 2) + pow(z, 2));
	}
	MatrixF Vector3D::GetRowVector() const
	{
		MatrixF retMat(3, 1);
		retMat.SetNum(0, 0, x);
		retMat.SetNum(1, 0, y);
		retMat.SetNum(2, 0, z);

		return retMat;
	}
	MatrixF Vector3D::GetColVector() const
	{
		MatrixF retMat(1, 3);
		retMat.SetNum(0, 0, x);
		retMat.SetNum(0, 1, y);
		retMat.SetNum(0, 2, z);

		return retMat;
	}
	Float3 Vector3D::TransformPoint(Float3 point)
	{
		point.x += x;
		point.y += y;
		point.z += z;

		return point;
	}

	Vector3D GetRelativeVec(const Vector3D& vec1, const Vector3D& vec2)
	{
		Vector3D retVec = Vector3D();
		retVec.x = vec2.x - vec1.x;
		retVec.y = vec2.y - vec1.y;
		retVec.z = vec2.z - (vec1.z - 1);

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
		return v1.x * v2.x + v1.y * v2.y;
	}
	Vector3D VectorCrossProduct(const Vector2D& v1, const Vector2D& v2)
	{
		MatrixF vectors = TransformF2x2(v1, v2);
		double zLength = MatrixGetDet(vectors);
		Vector3D cProduct = Vector3D(Float3(0, 0, zLength));

		return cProduct;
	}
	Vector3D VectorCrossProduct(const Vector3D& v1, const Vector3D& v2)
	{
		return Vector3D{
			v1.y * v2.z - v1.z * v2.y,
			v1.z * v2.x - v1.x * v2.z,
			v1.x * v2.y - v1.y * v2.x
		};
	}
	Vector2D VectorTripleProduct(const Vector2D& v1, const Vector2D& v2)
	{
		Vector3D tpVec = VectorCrossProduct(VectorCrossProduct(Vector3D(v1), Vector3D(v2)), Vector3D(v1));
		return Vector2D(tpVec.x, tpVec.y);
	}
	Vector3D VectorTripleProduct(const Vector3D& v1, const Vector3D& v2)
	{
		return VectorCrossProduct(VectorCrossProduct(Vector3D(v1), Vector3D(v2)), Vector3D(v1));
	}
	double VectorGetDeterminant(const Vector2D& v1, const Vector2D& v2)
	{
		MatrixF vectors = TransformF2x2(v1, v2);
		return MatrixGetDet(vectors);
	}
	void PrintProperties(const Vector2D& vector)
	{
		printf("X: %.3f, Y: %.3f\n", vector.x, vector.y);
	}
	void PrintProperties(const Vector3D& vector)
	{
		printf("X: %.3f, Y: %.3f, Z: %.3f\n", vector.x, vector.y, vector.z);
	}
}
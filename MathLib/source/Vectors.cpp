#include "../include/MathLib/MathLib.hpp"
#include <stdexcept>
#include <math.h>

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
	Vector2D::Vector2D(const float x, const float y)
		: Float2(x, y)
	{

	}
	Vector2D::Vector2D(const Line2D& line)
		: Float2(line.p2.x - line.p1.x, line.p2.y - line.p1.y)
	{
		
	}
	Vector2D::Vector2D(const Float2& p1, const Float2& p2)
		: Float2(p2.x - p1.x, p2.y - p1.y)
	{

	}
	const Vector2D& Vector2D::Transform(const Mat2& transformMat)
	{
		*this = *this * transformMat;
		return *this;
	}
	const Vector2D& Vector2D::Scale(const float scale)
	{
		x *= scale;
		y *= scale;
		return *this;
	}
	const Vector2D& Vector2D::Scale(const float scaleX, const float scaleY)
	{
		x *= scaleX;
		y *= scaleY;
		return *this;
	}
	const Vector2D& Vector2D::SetScale(const float scale)
	{
		float factor = sqrtf(Pow(x, 2) + Pow(y, 2));
		if (factor != 0)
			*this = *this * (float) (1.0 / factor);
		return *this;
	}
	const Vector2D& Vector2D::Rotate(float angle)
	{
		Mat2 rotMat = RotationMatrix2D(angle);

		*this = rotMat * *this;
		return *this;
	}
	const Vector2D& Vector2D::Rot90R()
	{
		float tmp = x;
		x = y;
		y = -tmp;
		return *this;
	}
	const Vector2D& Vector2D::Rot90L()
	{
		float tmp = x;
		x = -y;
		y = tmp;
		return *this;
	}
	const Vector2D& Vector2D::Normalize()
	{
		return SetScale(1);
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
		return Vector2D(x - other.x, y - other.y);
	}
	float Vector2D::operator*(const Vector2D& other) const
	{
		return x * other.x + y * other.y;
	}
	Vector2D Vector2D::operator*(const float number) const
	{
		return Vector2D(
			x * number,
			y * number
		);
	}
	Vector2D Vector2D::operator*(const Mat2& matrix) const
	{
		Mat2x1 result = GetRowVector() * matrix;
		return Vector2D(result.GetVal(0), result.GetVal(1));
	}
	void Vector2D::operator*=(const float number)
	{
		x *= number;
		y *= number;
	}
	Vector2D Vector2D::operator-() const
	{
		return Vector2D(-x, -y);
	}
	float Vector2D::GetAngle() const
	{
		float hyp = this->GetLen();
		float angle = fmod(180.001 + Rad2Deg(acos(y / hyp)), 180.001);

		if (x < 0)
			angle = -angle;

		return angle;
	}
	float Vector2D::GetLen() const
	{
		return sqrt(Pow(x, 2) + Pow(y, 2));
	}
	Mat2x1 Vector2D::GetRowVector() const
	{
		Mat2x1 retMat;
		retMat.SetVal(0, x);
		retMat.SetVal(1, y);

		return retMat;
	}
	Mat1x2 Vector2D::GetColVector() const
	{
		Mat1x2 retMat;
		retMat.SetVal(0, x);
		retMat.SetVal(1, y);

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
	Vector3D::Vector3D(const float x, const float y, const float z)
		: Float3(x, y, z)
	{

	}
	Vector3D::Vector3D(const Vector2D& vec, const float z)
		: Float3(vec.x, vec.y, z)
	{

	}
	Vector3D::Vector3D(const Float3& p1, const Float3& p2)
		: Float3(p2 - p1)
	{

	}
	const Vector3D& Vector3D::Transform(const Mat3& mat)
	{
		Mat3x1 result = GetRowVector() * mat;
		x = result.GetVal(0);
		y = result.GetVal(1);
		z = result.GetVal(2);
		return *this;
	}
	const Vector3D& Vector3D::Scale(const float& s)
	{
		x *= s;
		y *= s;
		z *= s;
		return *this;
	}
	const Vector3D& Vector3D::Scale(const float& sX, const float& sY, const float& sZ)
	{
		x *= sX;
		y *= sY;
		z *= sZ;
		return *this;
	}
	const Vector3D& Vector3D::SetScale(const float& len)
	{
		float vecLen = GetLen();

		Scale(vecLen * (len / vecLen));
		return *this;
	}
	const Vector3D& Vector3D::Rotate(const float& angle, const unsigned int& axis)
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
	const Vector3D& Vector3D::Normalize()
	{
		SetScale(1);
		return *this;
	}
	/*Vector3D Vector3D::Rotate(float angle, Vector3D *normal)
	{

	}*/
	Vector3D Vector3D::operator+(const Vector3D& other) const
	{
		return Vector3D(
			x + other.x,
			y + other.y,
			z + other.z
		);
	}
	Vector3D Vector3D::operator-(const Vector3D& other) const
	{
		return Vector3D(
			x - other.x,
			y - other.y,
			z - other.z
		);
	}
	float Vector3D::operator*(const Vector3D& other) const
	{
		return x * other.x + y * other.y + z * other.z;
	}
	Vector3D Vector3D::operator*(const float number) const
	{
		return Vector3D(x * number, y * number, z * number);
	}
	Vector3D Vector3D::operator*(const Mat3& matrix) const
	{
		Mat3x1 result = GetRowVector() * matrix;
		return Vector3D(result.GetVal(0), result.GetVal(1), result.GetVal(2));
	}
	void Vector3D::operator*=(const float number)
	{
		x *= number;
		y *= number;
		z *= number;
	}
	Vector3D Vector3D::operator-() const
	{
		return Vector3D(-x, -y, -z);
	}
	float Vector3D::GetAngle(const int axis) const
	{
		Vector2D calcVec = Vector2D();
		float angle;

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
	/*float Vector3D::GetAngle(Vector3D *normal)
	{

	}*/
	float Vector3D::GetLen() const
	{
		return sqrt(Pow(sqrt(Pow(x, 2) + Pow(y, 2)), 2) + Pow(z, 2));
	}
	Mat3x1 Vector3D::GetRowVector() const
	{
		Mat3x1 retMat;
		retMat.SetVal(0, x);
		retMat.SetVal(1, y);
		retMat.SetVal(2, z);

		return retMat;
	}
	Mat1x3 Vector3D::GetColVector() const
	{
		Mat1x3 retMat;
		retMat.SetVal(0, x);
		retMat.SetVal(1, y);
		retMat.SetVal(2, z);

		return retMat;
	}
	Float3 Vector3D::TransformPoint(Float3 point)
	{
		point.x += x;
		point.y += y;
		point.z += z;

		return point;
	}

	Vector2D operator*(const Mat2& matrix, const Vector2D& vector)
	{
		Mat1x2 colVec = vector.GetColVector();
		colVec = matrix * colVec;
		return Vector2D(colVec.GetVal(0), colVec.GetVal(1));
	}
	Vector3D operator*(const Mat3& matrix, const Vector3D& vector)
	{
		Mat1x3 colVec = vector.GetColVector();
		colVec = matrix * colVec;
		return Vector3D(colVec.GetVal(0), colVec.GetVal(1), colVec.GetVal(2));
	}
	Vector3D GetRelativeVec(const Vector3D& vec1, const Vector3D& vec2)
	{
		Vector3D retVec = Vector3D();
		retVec.x = vec2.x - vec1.x;
		retVec.y = vec2.y - vec1.y;
		retVec.z = vec2.z - (vec1.z - 1);

		return retVec;
	}
	float VectorGetAngleDifference(const Vector2D& v1, const Vector2D& v2)
	{
		float a1 = v1.GetAngle();
		float a2 = v2.GetAngle();

		float angle = (float) fmax(a1, a2) - (float) fmin(a1, a2);
		if (angle > 180)
			angle = 360 - angle;

		return angle;
	}
	float VectorDotProduct(const Vector2D& v1, const Vector2D& v2)
	{
		return v1.x * v2.x + v1.y * v2.y;
	}
	Vector3D VectorCrossProduct(const Vector2D& v1, const Vector2D& v2)
	{
		Mat2 vectors = Vec2Mat(v1, v2);
		float zLength = Determinant(vectors);
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
	float VectorGetDeterminant(const Vector2D& v1, const Vector2D& v2)
	{
		Mat2 vectors = Vec2Mat(v1, v2);
		return Determinant(vectors);
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
#include "../include/MathLib/MathLib.hpp"
#include <stdio.h>

namespace MathLib
{
	// Quaternion class

	Quaternion::Quaternion()
		: real(1), i(0), j(0), k(0)
	{

	}
	Quaternion::Quaternion(const float real, const float i, const float j, const float k)
		: real(real), i(i), j(j), k(k)
	{

	}
	Quaternion::Quaternion(const Float3& point)
		: real(0), i(point.x), j(point.y), k(point.z)
	{

	}
	Quaternion Quaternion::operator+(const Quaternion& other) const
	{
		return Quaternion(real + other.real, i + other.i, j + other.j, k + other.k);
	}
	Quaternion Quaternion::operator*(const Quaternion& other) const
	{
		Quaternion quat{};

		quat.real = real * other.real - i * other.i - j * other.j - k * other.k;
		quat.i = real * other.i + i * other.real + j * other.k - k * other.j;
		quat.j = real * other.j - i * other.k + j * other.real + k * other.i;
		quat.k = real * other.k + i * other.j - j * other.i + k * other.real;

		return quat;
	}
	void Quaternion::operator+=(const Quaternion& other)
	{
		real += other.real;
		i += other.i;
		j += other.j;
		k += other.k;
	}
	void Quaternion::operator*=(const Quaternion& other)
	{
		real = real * other.real - i * other.i - j * other.j - k * other.k;
		i = real * other.i + i * other.real + j * other.k - k * other.j;
		j = real * other.j - i * other.k + j * other.real + k * other.i;
		k = real * other.k + i * other.j - j * other.i + k * other.real;
	}
	Quaternion Quaternion::GetInverse() const
	{
		return Quaternion(real, -i, -j, -k);
	}
	Float3 Quaternion::GetPoint() const
	{
		return Float3(i, j, k);
	}
	Quaternion Quaternion::RotateQuaternion(const Quaternion& quat) const
	{
		return Quaternion(*this * quat * this->GetInverse());
	}
	Float3 Quaternion::RotatePoint(const Float3& point) const
	{
		Quaternion pQuat(point);
		return (*this * pQuat * this->GetInverse()).GetPoint();
	}

	Quaternion QuaternionRotation(const float angle, const float iAxis, const float jAxis, const float kAxis)
	{
		Quaternion retQuat(cos(Deg2Rad(angle)), 0, 0, 0);
		float sinAngle = sin(Deg2Rad(angle));
		retQuat.i = iAxis * sinAngle;
		retQuat.j = jAxis * sinAngle;
		retQuat.k = kAxis * sinAngle;
		return retQuat;
	}
	Quaternion QuaternionRotation(const float angle, const Float3& axis)
	{
		Quaternion retQuat(cos(Deg2Rad(angle / 2)), 0, 0, 0);
		float sinAngle = sin(Deg2Rad(angle / 2));

		float value = sqrtf(Square(axis.x) + Square(axis.y) + Square(axis.z));

		retQuat.i = axis.x / value * sinAngle;
		retQuat.j = axis.y / value * sinAngle;
		retQuat.k = axis.z / value * sinAngle;
		return retQuat;
	}
	void PrintProperties(const Quaternion& quat)
	{
		printf("Real: %.3f, i: %.3f, j: %.3f, k: %.3f\n", quat.real, quat.i, quat.j, quat.k);
	}

	// Operations

	Mat4 MatrixRotate(const Mat4& mat, const Quaternion& quat)
	{
		MathLib::Vector3D xVec{mat[0][0], mat[0][1], mat[0][2]};
		MathLib::Vector3D yVec{mat[1][0], mat[1][1], mat[1][2]};
		MathLib::Vector3D zVec{mat[2][0], mat[2][1], mat[2][2]};

		// A faster algorithm than pure Quaternion math
		// see https://blog.molecular-matters.com/2013/05/24/a-faster-quaternion-vector-multiplication/
		MathLib::Vector3D t = 2.0f * VectorCrossProduct(MathLib::Vector3D{quat.i, quat.j, quat.k}, xVec);
		xVec = MathLib::Float3(xVec) + quat.real * t + VectorCrossProduct(MathLib::Vector3D{quat.i, quat.j, quat.k}, t);

		t = 2.0f * VectorCrossProduct(MathLib::Vector3D{quat.i, quat.j, quat.k}, yVec);
		yVec = MathLib::Float3(yVec) + quat.real * t + VectorCrossProduct(MathLib::Vector3D{quat.i, quat.j, quat.k}, t);

		t = 2.0f * VectorCrossProduct(MathLib::Vector3D{quat.i, quat.j, quat.k}, zVec);
		zVec = MathLib::Float3(zVec) + quat.real * t + VectorCrossProduct(MathLib::Vector3D{quat.i, quat.j, quat.k}, t);

		Mat4 retMat{mat};

		retMat[0][0] = xVec.x;
		retMat[0][1] = xVec.y;
		retMat[0][2] = xVec.z;
		retMat[1][0] = yVec.x;
		retMat[1][1] = yVec.y;
		retMat[1][2] = yVec.z;
		retMat[2][0] = zVec.x;
		retMat[2][1] = zVec.y;
		retMat[2][2] = zVec.z;

		return retMat;
	}
}
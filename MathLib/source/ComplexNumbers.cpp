#include "../include/MathLib/MathLib.hpp"
#include <stdio.h>

namespace MathLib
{
	// Quaternion class

	Quaternion::Quaternion()
		: real(1), i(0), j(0), k(0)
	{

	}
	Quaternion::Quaternion(const double real, const double i, const double j, const double k)
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

	Quaternion QuaternionRotation(const double angle, const double iAxis, const double jAxis, const double kAxis)
	{
		Quaternion retQuat(cos(Deg2Rad(angle)), 0, 0, 0);
		double sinAngle = sin(Deg2Rad(angle));
		retQuat.i = iAxis * sinAngle;
		retQuat.j = jAxis * sinAngle;
		retQuat.k = kAxis * sinAngle;
		return retQuat;
	}
	Quaternion QuaternionRotation(const double angle, const Float3& axis)
	{
		Quaternion retQuat(cos(Deg2Rad(angle)), 0, 0, 0);
		double sinAngle = sin(Deg2Rad(angle));
		retQuat.i = axis.x * sinAngle;
		retQuat.j = axis.y * sinAngle;
		retQuat.k = axis.z * sinAngle;
		return retQuat;
	}
	void PrintProperties(const Quaternion& quat)
	{
		printf("Real: %.3f, i: %.3f, j: %.3f, k: %.3f\n", quat.real, quat.i, quat.j, quat.k);
	}

	// Operations

	Mat4 MatrixRotate(const Mat4& mat, const Quaternion& quat)
	{
		Quaternion xQuat{0, mat[0][0], mat[0][1], mat[0][2]};
		Quaternion yQuat{0, mat[1][0], mat[1][1], mat[1][2]};
		Quaternion zQuat{0, mat[2][0], mat[2][1], mat[2][2]};

		xQuat = quat.RotateQuaternion(xQuat);
		yQuat = quat.RotateQuaternion(yQuat);
		zQuat = quat.RotateQuaternion(zQuat);

		Mat4 retMat{mat};

		retMat[0][0] = xQuat.i;
		retMat[0][1] = xQuat.j;
		retMat[0][2] = xQuat.k;
		retMat[1][0] = yQuat.i;
		retMat[1][1] = yQuat.j;
		retMat[1][2] = yQuat.k;
		retMat[2][0] = zQuat.i;
		retMat[2][1] = zQuat.j;
		retMat[2][2] = zQuat.k;

		return retMat;
	}
}
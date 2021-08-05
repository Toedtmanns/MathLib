#include "../include/Maths.h"

namespace MathLib
{
	namespace Complex
	{
		imaginaryNum::imaginaryNum()
			: num(0)
		{

		}
		imaginaryNum::imaginaryNum(double num)
			: num(num)
		{

		}
		double imaginaryNum::operator+(const imaginaryNum &other)
		{
			return num + other.num;
		}
		double imaginaryNum::operator-(const imaginaryNum &other)
		{
			return num - other.num;
		}
		double imaginaryNum::operator*(const imaginaryNum &other)
		{
			if (num == other.num)
				return -1;
			else
				return num * other.num;
		}
		double imaginaryNum::operator/(const imaginaryNum &other)
		{
			return num / other.num;
		}
		bool imaginaryNum::operator==(const imaginaryNum &other)
		{
			return num == other.num;
		}
		bool imaginaryNum::operator!=(const imaginaryNum &other)
		{
			return num != other.num;
		}
		bool imaginaryNum::operator>=(const imaginaryNum &other)
		{
			return num >= other.num;
		}
		bool imaginaryNum::operator<=(const imaginaryNum &other)
		{
			return num <= other.num;
		}
		bool imaginaryNum::operator>(const imaginaryNum &other)
		{
			return num > other.num;
		}
		bool imaginaryNum::operator<(const imaginaryNum &other)
		{
			return num < other.num;
		}
		imaginaryNum& imaginaryNum::operator++()
		{
			++num;
			return *this;
		}
		imaginaryNum& imaginaryNum::operator--()
		{
			--num;
			return *this;
		}
		double imaginaryNum::operator++(int)
		{
			return num++;
		}
		double imaginaryNum::operator--(int)
		{
			return num--;
		}

		Quaternion::Quaternion()
			: real(1), i(0), j(0), k(0)
		{

		}
		Quaternion::Quaternion(double real, imaginaryNum i, imaginaryNum j, imaginaryNum k)
			: real(real), i(i), j(j), k(k)
		{

		}
		Quaternion Quaternion::operator+(const Quaternion &other)
		{
			return Quaternion(real + other.real, i + other.i, j + other.j, k + other.k);
		}
		Quaternion Quaternion::operator*(const Quaternion &other)
		{
			return Quaternion(
				real * other.real - i * other.i - j * other.j - k - other.k,
				real * other.i + i * other.real + j * other.k - k * other.j,
				real * other.j - i * other.k + j * other.real + k * other.i,
				real * other.k + i * other.j - j * other.i + k * other.real
			);
		}
	}
}
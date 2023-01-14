#pragma once
#include "mlPrimitives.hpp"
#include "mlMatrices.hpp"

namespace MathLib
{
	// Vector math

	class Vector2D : public Float2
	{
	public:
		constexpr Vector2D()
			: Float2(0, 0)
		{

		}
		constexpr Vector2D(const Float2& dir)
			: Float2(dir)
		{

		}
		constexpr explicit Vector2D(const float x, const float y)
			: Float2(x, y)
		{

		}
		constexpr Vector2D(const Line2D& line)
			: Float2(line.p2.x - line.p1.x, line.p2.y - line.p1.y)
		{

		}
		constexpr Vector2D(const Float2& p1, const Float2& p2)
			: Float2(p2.x - p1.x, p2.y - p1.y)
		{

		}
		
		constexpr Vector2D operator+(const Vector2D& other) const
		{
			return Vector2D(Float2(x + other.x, y + other.y));
		}
		constexpr Vector2D operator-(const Vector2D& other) const
		{
			return Vector2D(x - other.x, y - other.y);
		}
		constexpr float operator*(const Vector2D& other) const
		{
			return x * other.x + y * other.y;
		}
		constexpr Vector2D operator*(const float number) const
		{
			return Vector2D(
				x * number,
				y * number
			);
		}
		constexpr Vector2D operator*(const Mat2& matrix) const
		{
			Mat2x1 result = GetRowVector() * matrix;
			return Vector2D(result.GetVal(0), result.GetVal(1));
		}
		constexpr void operator*=(const float number)
		{
			x *= number;
			y *= number;
		}
		constexpr Vector2D operator-() const
		{
			return Vector2D(-x, -y);
		}

		inline const Vector2D& Transform(const Mat2& transformMat);
		constexpr const Vector2D& Scale(const float scale)
		{
			x *= scale;
			y *= scale;
			return *this;
		}
		constexpr const Vector2D& Scale(const float scaleX, const float scaleY)
		{
			x *= scaleX;
			y *= scaleY;
			return *this;
		}
		const Vector2D& SetScale(const float scale)
		{
			float factor = sqrtf(Pow(x, 2) + Pow(y, 2));
			if (factor != 0)
				*this = *this * (float) (1.0 / factor);
			return *this;
		}
		constexpr const Vector2D& Rotate(float angle)
		{
			Mat2 rotMat = RotationMatrix2D(angle);

			*this = *this * rotMat;
			return *this;
		}
		constexpr const Vector2D& Rot90R()
		{
			float tmp = x;
			x = y;
			y = -tmp;
			return *this;
		}
		constexpr const Vector2D& Rot90L()
		{
			float tmp = x;
			x = -y;
			y = tmp;
			return *this;
		}
		const Vector2D& Normalize()
		{
			return SetScale(1);
		}
		constexpr Vector2D TripleProduct(const Vector2D& other) const;

		inline float GetAngle() const
		{
			float hyp = this->GetLen();
			float angle = fmod(180.001 + Rad2Deg(acos(y / hyp)), 180.001);

			if (x < 0)
				angle = -angle;

			return angle;
		}
		inline float GetLen() const
		{
			return sqrt(Pow(x, 2) + Pow(y, 2));
		}
		constexpr Mat2x1 GetRowVector() const
		{
			Mat2x1 retMat;
			retMat.SetVal(0, x);
			retMat.SetVal(1, y);

			return retMat;
		}
		constexpr Mat1x2 GetColVector() const
		{
			Mat1x2 retMat;
			retMat.SetVal(0, x);
			retMat.SetVal(1, y);

			return retMat;
		}
		constexpr Float2 TransformPoint(Float2 point) const
		{
			point.x += x;
			point.y += y;

			return point;
		}
	};

	class Vector3D : public Float3
	{
	public:
		constexpr Vector3D()
			: Float3(0, 0, 0)
		{

		}
		constexpr Vector3D(const Float3& dir)
			: Float3(dir)
		{

		}
		constexpr explicit Vector3D(const float x, const float y, const float z)
			: Float3(x, y, z)
		{

		}
		constexpr Vector3D(const Line3D& line)
			: Float3(line.p2 - line.p1)
		{

		}
		constexpr Vector3D(const Vector2D& vec, const float z = 0)
			: Float3(vec.x, vec.y, z)
		{

		}
		constexpr Vector3D(const Float3& p1, const Float3& p2)
			: Float3(p2 - p1)
		{

		}
		constexpr const Vector3D& Transform(const Mat3& mat)
		{
			Mat3x1 result = GetRowVector() * mat;
			x = result.GetVal(0);
			y = result.GetVal(1);
			z = result.GetVal(2);
			return *this;
		}
		constexpr const Vector3D& Scale(const float& s)
		{
			x *= s;
			y *= s;
			z *= s;
			return *this;
		}
		constexpr const Vector3D& Scale(const float& sX, const float& sY, const float& sZ)
		{
			x *= sX;
			y *= sY;
			z *= sZ;
			return *this;
		}
		const Vector3D& SetScale(const float& len)
		{
			float vecLen = GetLen();

			Scale(len / vecLen);
			return *this;
		}
		constexpr const Vector3D& Rotate(const float& angle, const unsigned int& axis)
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
		inline const Vector3D& Normalize()
		{
			SetScale(1);
			return *this;
		}
		constexpr Vector3D operator+(const Vector3D& other) const
		{
			return Vector3D(
				x + other.x,
				y + other.y,
				z + other.z
			);
		}
		constexpr Vector3D operator-(const Vector3D& other) const
		{
			return Vector3D(
				x - other.x,
				y - other.y,
				z - other.z
			);
		}
		constexpr float operator*(const Vector3D& other) const
		{
			return x * other.x + y * other.y + z * other.z;
		}
		constexpr Vector3D operator*(const float number) const
		{
			return Vector3D(x * number, y * number, z * number);
		}
		constexpr Vector3D operator*(const Mat3& matrix) const
		{
			Mat3x1 result = GetRowVector() * matrix;
			return Vector3D(result.GetVal(0), result.GetVal(1), result.GetVal(2));
		}
		constexpr void operator*=(const float number)
		{
			x *= number;
			y *= number;
			z *= number;
		}
		constexpr Vector3D operator-() const
		{
			return Vector3D(-x, -y, -z);
		}
		float GetAngle(const int axis) const
		{
			Vector2D calcVec = Vector2D();
			float angle = 0.0f;

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
			}

			return angle;
		}
		float GetLen() const
		{
			return sqrt(Pow(x, 2) + Pow(y, 2) + Pow(z, 2));
		}
		constexpr Mat3x1 GetRowVector() const
		{
			Mat3x1 retMat;
			retMat.SetVal(0, x);
			retMat.SetVal(1, y);
			retMat.SetVal(2, z);

			return retMat;
		}
		constexpr Mat1x3 GetColVector() const
		{
			Mat1x3 retMat;
			retMat.SetVal(0, x);
			retMat.SetVal(1, y);
			retMat.SetVal(2, z);

			return retMat;
		}
		constexpr bool IsNull() const
		{
			return x == 0.0 && y == 0.0 && z == 0.0;
		}
		constexpr Float3 TransformPoint(Float3 point)
		{
			point.x += x;
			point.y += y;
			point.z += z;

			return point;
		}
	};

	// Operator overloads

	inline Vector2D operator*(const Mat2& matrix, const Vector2D& vector)
	{
		Mat1x2 colVec = vector.GetColVector();
		colVec = matrix * colVec;
		return Vector2D(colVec.GetVal(0), colVec.GetVal(1));
	}
	inline Vector3D operator*(const Mat3& matrix, const Vector3D& vector)
	{
		Mat1x3 colVec = vector.GetColVector();
		colVec = matrix * colVec;
		return Vector3D(colVec.GetVal(0), colVec.GetVal(1), colVec.GetVal(2));
	}

	// Other utility functions

	constexpr Vector3D GetRelativeVec(const Vector3D& vec1, const Vector3D& vec2)
	{
		Vector3D retVec = Vector3D();
		retVec.x = vec2.x - vec1.x;
		retVec.y = vec2.y - vec1.y;
		retVec.z = vec2.z - (vec1.z - 1);

		return retVec;
	}
	inline float VectorGetAngleDifference(const Vector2D& v1, const Vector2D& v2)
	{
		float a1 = v1.GetAngle();
		float a2 = v2.GetAngle();

		float angle = (float) Max(a1, a2) - (float) Min(a1, a2);
		if (angle > 180)
			angle = 360 - angle;

		return angle;
	}
	constexpr float VectorDotProduct(const Vector2D& v1, const Vector2D& v2)
	{
		return v1.x * v2.x + v1.y * v2.y;
	}
	constexpr float VectorDotProduct(const Vector3D& v1, const Vector3D& v2)
	{
		return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
	}
	constexpr Vector3D VectorCrossProduct(const Vector2D& v1, const Vector2D& v2)
	{
		Mat2 vectors = Vec2Mat(v1, v2);
		float zLength = Determinant(vectors);
		Vector3D cProduct = Vector3D(Float3(0, 0, zLength));

		return cProduct;
	}
	constexpr Vector3D VectorCrossProduct(const Vector3D& v1, const Vector3D& v2)
	{
		return Vector3D{
			v1.y * v2.z - v1.z * v2.y,
			v1.z * v2.x - v1.x * v2.z,
			v1.x * v2.y - v1.y * v2.x
		};
	}
	constexpr Vector2D VectorTripleProduct(const Vector2D& v1, const Vector2D& v2)
	{
		Vector3D tpVec = VectorCrossProduct(VectorCrossProduct(Vector3D(v1), Vector3D(v2)), Vector3D(v1));
		return Vector2D(tpVec.x, tpVec.y);
	}
	constexpr Vector3D VectorTripleProduct(const Vector3D& v1, const Vector3D& v2)
	{
		return VectorCrossProduct(VectorCrossProduct(Vector3D(v1), Vector3D(v2)), Vector3D(v1));
	}
	constexpr float VectorGetDeterminant(const Vector2D& v1, const Vector2D& v2)
	{
		Mat2 vectors = Vec2Mat(v1, v2);
		return Determinant(vectors);
	}
	inline Mat3 Vec2Mat(const Vector3D& vector)
	{
		Vector3D perp1 = VectorCrossProduct(vector, Vector3D(-vector.z, vector.x, vector.y)).Normalize();
		Vector3D perp2 = VectorCrossProduct(vector, perp1).Normalize();
		float values[9] = {
			perp1.x, perp1.y, perp1.z,
			perp2.x, perp2.y, perp2.z,
			vector.x, vector.y, vector.z
		};

		Mat3 mat;
		mat.SetMatrix(values);
		return mat;
	}

	// Some late definitions of methods of the vector classes

	inline const Vector2D& Vector2D::Transform(const Mat2& transformMat)
	{
		*this = transformMat * *this;
		return *this;
	}
	constexpr Vector2D Vector2D::TripleProduct(const Vector2D& other) const
	{
		Vector3D tpVec = VectorCrossProduct(VectorCrossProduct(Vector3D(*this), Vector3D(other)), Vector3D(*this));
		return Vector2D(tpVec.x, tpVec.y);
	}
}
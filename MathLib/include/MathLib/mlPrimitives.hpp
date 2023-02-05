#pragma once
#include "mlCommon.hpp"

namespace MathLib
{
	// Template primitive
	template <typename Type>
	class Primitive2
	{
	public:
		union
		{
			struct
			{
				Type x;
				Type y;
			};
			Type content[2];
		};

		constexpr Primitive2()
			: x(0), y(0)
		{

		}
		constexpr Primitive2(const Type x, const Type y)
			: x(x), y(y)
		{

		}
		constexpr Primitive2(const Primitive2<Type>& other)
		{
			x = (Type) other.x;
			y = (Type) other.y;
		}
		constexpr void operator=(const Primitive2<Type>& other)
		{
			x = (Type) other.x;
			y = (Type) other.y;
		}
		template <typename OtherType>
		constexpr Primitive2(const Primitive2<OtherType>& other)
		{
			y = (Type) other.y;
			x = (Type) other.x;
		}
		template <typename OtherType>
		constexpr void operator=(const Primitive2<OtherType>& other)
		{
			x = (Type) other.x;
			y = (Type) other.y;
		}
		constexpr Type operator[](const size_t index)
		{
			return content[index];
		}
		constexpr const Type operator[](const size_t index) const
		{
			return content[index];
		}
		constexpr Primitive2<Type> operator-() const
		{
			return Primitive2<Type>(-x, -y);
		}
	};

	template <typename Type>
	class Primitive3
	{
	public:
		union
		{
			struct
			{
				Type x;
				Type y;
				Type z;
			};
			struct
			{
				Type r;
				Type g;
				Type b;
			};
			Type content[3];
		};

		constexpr Primitive3()
			: x(0), y(0), z(0)
		{

		}
		constexpr Primitive3(const Type x, const Type y, const Type z)
			: x(x), y(y), z(z)
		{

		}
		constexpr Primitive3(const Primitive3<Type>& other)
		{
			x = (Type) other.x;
			y = (Type) other.y;
			z = (Type) other.z;
		}
		constexpr void operator=(const Primitive3<Type>& other)
		{
			x = (Type) other.x;
			y = (Type) other.y;
			z = (Type) other.z;
		}
		template <typename OtherType>
		constexpr Primitive3(const Primitive3<OtherType>& other)
		{
			x = (Type) round(other.x);
			y = (Type) round(other.y);
			z = (Type) round(other.z);
		}
		template <typename OtherType>
		constexpr void operator=(const Primitive3<OtherType>& other)
		{
			x = (Type) round(other.x);
			y = (Type) round(other.y);
			z = (Type) round(other.z);
		}
		constexpr Type operator[](const size_t index)
		{
			return content[index];
		}
		constexpr const Type operator[](const size_t index) const
		{
			return content[index];
		}
		constexpr Primitive3<Type> operator-() const
		{
			return Primitive3<Type>(-x, -y, -z);
		}
	};

	template <typename Type>
	class Primitive4
	{
	public:
		union
		{
			struct
			{
				Type x;
				Type y;
				Type z;
				Type w;
			};
			struct
			{
				Type r;
				Type g;
				Type b;
				Type a;
			};
			Type content[4];
		};

		constexpr Primitive4()
			: x(0), y(0), z(0), w(0)
		{

		}
		constexpr Primitive4(const Type x, const Type y, const Type z, const Type w)
			: x(x), y(y), z(z), w(w)
		{

		}
		constexpr Primitive4(const Primitive4<Type>& other)
		{
			x = (Type) other.x;
			y = (Type) other.y;
			z = (Type) other.z;
			w = (Type) other.w;
		}
		constexpr void operator=(const Primitive4<Type>& other)
		{
			x = (Type) other.x;
			y = (Type) other.y;
			z = (Type) other.z;
			w = (Type) other.w;
		}
		template <typename OtherType>
		constexpr Primitive4(const Primitive4<OtherType>& other)
		{
			x = (Type) round(other.x);
			y = (Type) round(other.y);
			z = (Type) round(other.z);
			w = (Type) round(other.w);
		}
		template <typename OtherType>
		constexpr void operator=(const Primitive4<OtherType>& other)
		{
			x = (Type) round(other.x);
			y = (Type) round(other.y);
			z = (Type) round(other.z);
			w = (Type) round(other.w);
		}
		constexpr Type operator[](const size_t index)
		{
			return content[index];
		}
		constexpr const Type operator[](const size_t index) const
		{
			return content[index];
		}
		constexpr Primitive4<Type> operator-() const
		{
			return Primitive4<Type>(-x, -y, -z, -w);
		}
	};

	// Primitive types

	typedef Primitive2<int> Int2;
	typedef Primitive3<int> Int3;
	typedef Primitive4<int> Int4;

	typedef Primitive2<float> Float2;
	typedef Primitive3<float> Float3;
	typedef Primitive4<float> Float4;


	// Operator overloads

	template <typename Type>
	constexpr bool operator==(const Primitive2<Type>& f1, const Primitive2<Type>& f2)
	{
		return f1.x == f2.x && f1.y == f2.y;
	}
	template <typename Type>
	constexpr bool operator!=(const Primitive2<Type>& f1, const Primitive2<Type>& f2)
	{
		return f1.x != f2.x || f1.y != f2.y;
	}
	template <typename Type>
	constexpr bool operator==(const Primitive3<Type>& f1, const Primitive3<Type>& f2)
	{
		return f1.x == f2.x && f1.y == f2.y && f1.z == f2.z;
	}
	template <typename Type>
	constexpr bool operator!=(const Primitive3<Type>& f1, const Primitive3<Type>& f2)
	{
		return f1.x != f2.x || f1.y != f2.y || f1.z != f2.z;
	}
	template <typename Type>
	constexpr bool operator==(const Primitive4<Type>& f1, const Primitive4<Type>& f2)
	{
		return f1.x == f2.x && f1.y == f2.y && f1.z == f2.z && f1.w == f2.w;
	}
	template <typename Type>
	constexpr bool operator!=(const Primitive4<Type> & f1, const Primitive4<Type>& f2)
	{
		return f1.x != f2.x || f1.y != f2.y || f1.z != f2.z || f1.w != f2.w;
	}
	template <typename Type>
	constexpr Primitive2<Type> operator+(const Primitive2<Type>& f1, const Primitive2<Type>& f2)
	{
		return Primitive2<Type>(f1.x + f2.x, f1.y + f2.y);
	}
	template <typename Type>
	constexpr Primitive3<Type> operator+(const Primitive3<Type>& f1, const Primitive3<Type>& f2)
	{
		return Primitive3<Type>(f1.x + f2.x, f1.y + f2.y, f1.z + f2.z);
	}
	template <typename Type>
	constexpr Primitive4<Type> operator+(const Primitive4<Type>& f1, const Primitive4<Type>& f2)
	{
		return Primitive4<Type>(f1.x + f2.x, f1.y + f2.y, f1.z + f2.z, f1.w + f2.w);
	}
	template <typename Type>
	constexpr Primitive2<Type> operator-(const Primitive2<Type>& f1, const Primitive2<Type>& f2)
	{
		return Primitive2<Type>(f1.x - f2.x, f1.y - f2.y);
	}
	template <typename Type>
	constexpr Primitive3<Type> operator-(const Primitive3<Type>& f1, const Primitive3<Type>& f2)
	{
		return Primitive3<Type>(f1.x - f2.x, f1.y - f2.y, f1.z - f2.z);
	}
	template <typename Type>
	constexpr Primitive4<Type> operator-(const Primitive4<Type>& f1, const Primitive4<Type>& f2)
	{
		return Primitive4<Type>(f1.x - f2.x, f1.y - f2.y, f1.z - f2.z, f1.w - f2.w);
	}
	template <typename Type>
	constexpr void operator+=(Primitive2<Type>& f1, const Primitive2<Type>& f2)
	{
		f1 = Primitive2<Type>(f1.x + f2.x, f1.y + f2.y);
	}
	template <typename Type>
	constexpr void operator+=(Primitive3<Type>& f1, const Primitive3<Type>& f2)
	{
		f1 = Primitive3<Type>(f1.x + f2.x, f1.y + f2.y, f1.z + f2.z);
	}
	template <typename Type>
	constexpr void operator+=(Primitive4<Type>& f1, const Primitive4<Type>& f2)
	{
		f1 = Primitive4<Type>(f1.x + f2.x, f1.y + f2.y, f1.z + f2.z, f1.w + f2.w);
	}
	template <typename Type>
	constexpr void operator-=(Primitive2<Type>& f1, const Primitive2<Type>& f2)
	{
		f1 = Primitive2<Type>(f1.x - f2.x, f1.y - f2.y);
	}
	template <typename Type>
	constexpr void operator-=(Primitive3<Type>& f1, const Primitive3<Type>& f2)
	{
		f1 = Primitive3<Type>(f1.x - f2.x, f1.y - f2.y, f1.z - f2.z);
	}
	template <typename Type>
	constexpr void operator-=(Primitive4<Type>& f1, const Primitive4<Type>& f2)
	{
		f1 = Primitive4<Type>(f1.x - f2.x, f1.y - f2.y, f1.z - f2.z, f1.w - f2.w);
	}
	template <typename Type>
	constexpr Primitive2<Type> operator*(const Primitive2<Type>& f1, const Primitive2<Type>& f2)
	{
		return Primitive2<Type>(f1.x * f2.x, f1.y * f2.y);
	}
	template <typename Type>
	constexpr Primitive3<Type> operator*(const Primitive3<Type>& f1, const Primitive3<Type>& f2)
	{
		return Primitive3<Type>(f1.x * f2.x, f1.y * f2.y, f1.z * f2.z);
	}
	template <typename Type>
	constexpr Primitive4<Type> operator*(const Primitive4<Type>& f1, const Primitive4<Type>& f2)
	{
		return Primitive4<Type>(f1.x * f2.x, f1.y * f2.y, f1.z * f2.z, f1.w * f2.w);
	}
	template <typename Type>
	constexpr Primitive2<Type> operator*(const Primitive2<Type>& f1, const double f2)
	{
		return Primitive2<Type>(f1.x * f2, f1.y * f2);
	}
	template <typename Type>
	constexpr Primitive3<Type> operator*(const Primitive3<Type>& f1, const double f2)
	{
		return Primitive3<Type>(f1.x * f2, f1.y * f2, f1.z * f2);
	}
	template <typename Type>
	constexpr Primitive4<Type> operator*(const Primitive4<Type>& f1, const double f2)
	{
		return Primitive4<Type>(f1.x * f2, f1.y * f2, f1.z * f2, f1.w * f2);
	}
	template <typename Type>
	constexpr Primitive2<Type> operator*(const double f2, const Primitive2<Type>& f1)
	{
		return Primitive2<Type>(f1.x * f2, f1.y * f2);
	}
	template <typename Type>
	constexpr Primitive3<Type> operator*(const double f2, const Primitive3<Type>& f1)
	{
		return Primitive3<Type>(f1.x * f2, f1.y * f2, f1.z * f2);
	}
	template <typename Type>
	constexpr Primitive4<Type> operator*(const double f2, const Primitive4<Type>& f1)
	{
		return Primitive4<Type>(f1.x * f2, f1.y * f2, f1.z * f2, f1.w * f2);
	}

	template <typename Type>
	constexpr Primitive2<Type> operator*=(Primitive2<Type>& f1, const Primitive2<Type>& f2)
	{
		f1 = Primitive2<Type>(f1.x * f2.x, f1.y * f2.y);
	}
	template <typename Type>
	constexpr Primitive3<Type> operator*=(Primitive3<Type>& f1, const Primitive3<Type>& f2)
	{
		f1 = Primitive3<Type>(f1.x * f2.x, f1.y * f2.y, f1.z * f2.z);
	}
	template <typename Type>
	constexpr Primitive4<Type> operator*=(Primitive4<Type>& f1, const Primitive4<Type>& f2)
	{
		f1 = Primitive4<Type>(f1.x * f2.x, f1.y * f2.y, f1.z * f2.z, f1.w * f2.w);
	}
	template <typename Type>
	constexpr Primitive2<Type> operator*=(Primitive2<Type>& f1, const double f2)
	{
		f1 = Primitive2<Type>(f1.x * f2, f1.y * f2);
	}
	template <typename Type>
	constexpr Primitive3<Type> operator*=(Primitive3<Type>& f1, const double f2)
	{
		f1 = Primitive3<Type>(f1.x * f2, f1.y * f2, f1.z * f2);
	}
	template <typename Type>
	constexpr Primitive4<Type> operator*=(Primitive4<Type>& f1, const double f2)
	{
		f1 = Primitive4<Type>(f1.x * f2, f1.y * f2, f1.z * f2, f1.w * f2);
	}

	// Other utility functions

	template <typename Type>
	constexpr Primitive2<Type> Lerp(const Primitive2<Type>& start, const Primitive2<Type>& end, const float t)
	{
		return {
			start.x + t * (end.x - start.x),
			start.y + t * (end.y - start.y)
		};
	}
	template <typename Type>
	constexpr Primitive3<Type> Lerp(const Primitive3<Type>& start, const Primitive3<Type>& end, const float t)
	{
		return {
			start.x + t * (end.x - start.x),
			start.y + t * (end.y - start.y),
			start.z + t * (end.z - start.z)
		};
	}
	template <typename Type>
	constexpr Primitive4<Type> Lerp(const Primitive4<Type>& start, const Primitive4<Type>& end, const float t)
	{
		return {
			start.x + t * (end.x - start.x),
			start.y + t * (end.y - start.y),
			start.z + t * (end.z - start.z),
			start.w + t * (end.w - start.w)
		};
	}
	template <typename Type>
	float GetDistance(const Primitive2<Type>& p1, const Primitive2<Type>& p2)
	{
		return sqrt(Pow(p2.x - p1.x, 2) + Pow(p2.y - p1.y, 2));
	}
	template <typename Type>
	float GetDistance(const Primitive3<Type>& p1, const Primitive3<Type>& p2)
	{
		return sqrt(Pow(p2.x - p1.x, 2) + Pow(p2.y - p1.y, 2) + Pow(p2.z - p1.z, 2));
	}
	template <typename Type>
	float GetDistance(const Primitive4<Type>& p1, const Primitive4<Type>& p2)
	{
		return sqrt(Pow(p2.x - p1.x, 2) + Pow(p2.y - p1.y, 2) + Pow(p2.z - p1.z, 2) + Pow(p2.w - p1.w, 2));
	}

	// Lines

	struct Line2D
	{
		Float2 p1;
		Float2 p2;

		constexpr Line2D()
			: p1(), p2()
		{

		}
		constexpr Line2D(const float x1, const float y1, const float x2, const float y2)
			: p1(x1, y1), p2(x2, y2)
		{

		}
		constexpr Line2D(const Float2& p1, const Float2& p2)
			: p1(p1), p2(p2)
		{

		}
		constexpr Float2 GetInnerDistance() const
		{
			float xDist = p1.x - p2.x;
			float yDist = p1.y - p2.y;

			if (xDist < 0)
				xDist = -xDist;
			if (yDist < 0)
				yDist = -yDist;

			return Float2(xDist, yDist);
		}
		float GetLength() const
		{
			return sqrt(Pow(p2.x - p1.x, 2) * Pow(p2.y - p1.y, 2));
		}
		constexpr bool operator==(const Line2D& other)
		{
			return p1 == other.p1 && p2 == other.p2;
		}
		constexpr bool operator!=(const Line2D& other)
		{
			return !(p1 == other.p1 && p2 == other.p2);
		}
	};

	struct Line3D
	{
		Float3 p1;
		Float3 p2;

		constexpr Line3D()
			: p1(), p2()
		{

		}
		constexpr Line3D(const Float3& p1, const Float3& p2)
			: p1(p1), p2(p2)
		{

		}
		constexpr Float3 GetInnerDistance() const
		{
			float xDist = p1.x - p2.x;
			float yDist = p1.y - p2.y;
			float zDist = p1.z - p2.z;

			if (xDist < 0)
				xDist = -xDist;
			if (yDist < 0)
				yDist = -yDist;
			if (zDist < 0)
				zDist = -zDist;

			return Float3(xDist, yDist, zDist);
		}
		float GetLength() const
		{
			return GetDistance(p1, p2);
		}
		constexpr bool operator==(const Line3D& other)
		{
			return p1 == other.p1 && p2 == other.p2;
		}
		constexpr bool operator!=(const Line3D& other)
		{
			return !(p1 == other.p1 && p2 == other.p2);
		}
	};
}
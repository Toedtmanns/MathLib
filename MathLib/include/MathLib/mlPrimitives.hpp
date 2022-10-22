#pragma once
#include "mlCommon.hpp"

namespace MathLib
{
	// Template primitive
	template <typename Type>
	class EXPORT Primitive2
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
		Primitive2(const Primitive2<Type>& other)
		{
			x = (Type) other.x;
			y = (Type) other.y;
		}
		void operator=(const Primitive2<Type>& other)
		{
			x = (Type) other.x;
			y = (Type) other.y;
		}
		template <typename OtherType>
		Primitive2(const Primitive2<OtherType>& other)
		{
			x = (Type) round(other.x);
			y = (Type) round(other.y);
		}
		template <typename OtherType>
		void operator=(const Primitive2<OtherType>& other)
		{
			x = (Type) round(other.x);
			y = (Type) round(other.y);
		}
		constexpr double operator[](const size_t index)
		{
			return content[index];
		}
		constexpr const double operator[](const size_t index) const
		{
			return content[index];
		}
	};

	template <typename Type>
	class EXPORT Primitive3
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
		Primitive3(const Primitive3<Type>& other)
		{
			x = (Type) other.x;
			y = (Type) other.y;
			z = (Type) other.z;
		}
		void operator=(const Primitive3<Type>& other)
		{
			x = (Type) other.x;
			y = (Type) other.y;
			z = (Type) other.z;
		}
		template <typename OtherType>
		Primitive3(const Primitive3<OtherType>& other)
		{
			x = (Type) round(other.x);
			y = (Type) round(other.y);
			z = (Type) round(other.z);
		}
		template <typename OtherType>
		void operator=(const Primitive3<OtherType>& other)
		{
			x = (Type) round(other.x);
			y = (Type) round(other.y);
			z = (Type) round(other.z);
		}
		constexpr double operator[](const size_t index)
		{
			return content[index];
		}
		constexpr const double operator[](const size_t index) const
		{
			return content[index];
		}
	};

	template <typename Type>
	class EXPORT Primitive4
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
			: x(0), y(0), z(0)
		{

		}
		constexpr Primitive4(const Type x, const Type y, const Type z, const Type w)
			: x(x), y(y), z(z), w(w)
		{

		}
		Primitive4(const Primitive4<Type>& other)
		{
			x = (Type) other.x;
			y = (Type) other.y;
			z = (Type) other.z;
			w = (Type) other.w;
		}
		void operator=(const Primitive4<Type>& other)
		{
			x = (Type) other.x;
			y = (Type) other.y;
			z = (Type) other.z;
			w = (Type) other.w;
		}
		template <typename OtherType>
		Primitive4(const Primitive4<OtherType>& other)
		{
			x = (Type) round(other.x);
			y = (Type) round(other.y);
			z = (Type) round(other.z);
			w = (Type) round(other.w);
		}
		template <typename OtherType>
		void operator=(const Primitive4<OtherType>& other)
		{
			x = (Type) round(other.x);
			y = (Type) round(other.y);
			z = (Type) round(other.z);
			w = (Type) round(other.w);
		}
		constexpr double operator[](const size_t index)
		{
			return content[index];
		}
		constexpr const double operator[](const size_t index) const
		{
			return content[index];
		}
	};

	// Primitive types

	typedef Primitive2<int> Int2;
	typedef Primitive3<int> Int3;
	typedef Primitive4<int> Int4;

	typedef Primitive2<double> Float2;
	typedef Primitive3<double> Float3;
	typedef Primitive4<double> Float4;


	// Operator overloads

	template <typename Type>
	EXPORT constexpr bool operator==(const Primitive2<Type>& f1, const Primitive2<Type>& f2)
	{
		return f1.x == f2.x && f1.y == f2.y;
	}
	template <typename Type>
	EXPORT constexpr bool operator!=(const Primitive2<Type>& f1, const Primitive2<Type>& f2)
	{
		return f1.x != f2.x || f1.y != f2.y;
	}
	template <typename Type>
	EXPORT constexpr bool operator==(const Primitive3<Type>& f1, const Primitive3<Type>& f2)
	{
		return f1.x == f2.x && f1.y == f2.y && f1.z == f2.z;
	}
	template <typename Type>
	EXPORT constexpr bool operator!=(const Primitive3<Type>& f1, const Primitive3<Type>& f2)
	{
		return f1.x != f2.x || f1.y != f2.y || f1.z != f2.z;
	}
	template <typename Type>
	EXPORT constexpr bool operator==(const Primitive4<Type>& f1, const Primitive4<Type>& f2)
	{
		return f1.x == f2.x && f1.y == f2.y && f1.z == f2.z && f1.w == f2.w;
	}
	template <typename Type>
	EXPORT constexpr bool operator!=(const Primitive4<Type> & f1, const Primitive4<Type>& f2)
	{
		return f1.x != f2.x || f1.y != f2.y || f1.z != f2.z || f1.w != f2.w;
	}
	template <typename Type>
	EXPORT constexpr Primitive2<Type> operator+(const Primitive2<Type>& f1, const Primitive2<Type>& f2)
	{
		return Primitive2<Type>(f1.x + f2.x, f1.y + f2.y);
	}
	template <typename Type>
	EXPORT constexpr Primitive3<Type> operator+(const Primitive3<Type>& f1, const Primitive3<Type>& f2)
	{
		return Primitive3<Type>(f1.x + f2.x, f1.y + f2.y, f1.z + f2.z);
	}
	template <typename Type>
	EXPORT constexpr Primitive4<Type> operator+(const Primitive4<Type>& f1, const Primitive4<Type>& f2)
	{
		return Primitive4<Type>(f1.x + f2.x, f1.y + f2.y, f1.z + f2.z, f1.w + f2.w);
	}
	template <typename Type>
	EXPORT constexpr Primitive2<Type> operator-(const Primitive2<Type>& f1, const Primitive2<Type>& f2)
	{
		return Primitive2<Type>(f1.x - f2.x, f1.y - f2.y);
	}
	template <typename Type>
	EXPORT constexpr Primitive3<Type> operator-(const Primitive3<Type>& f1, const Primitive3<Type>& f2)
	{
		return Primitive3<Type>(f1.x - f2.x, f1.y - f2.y, f1.z - f2.z);
	}
	template <typename Type>
	EXPORT constexpr Primitive4<Type> operator-(const Primitive4<Type>& f1, const Primitive4<Type>& f2)
	{
		return Primitive4<Type>(f1.x - f2.x, f1.y - f2.y, f1.z - f2.z, f1.w - f2.w);
	}
	template <typename Type>
	EXPORT constexpr void operator+=(Primitive2<Type>& f1, const Primitive2<Type>& f2)
	{
		f1 = Primitive2<Type>(f1.x + f2.x, f1.y + f2.y);
	}
	template <typename Type>
	EXPORT constexpr void operator+=(Primitive3<Type>& f1, const Primitive3<Type>& f2)
	{
		f1 = Primitive3<Type>(f1.x + f2.x, f1.y + f2.y, f1.z + f2.z);
	}
	template <typename Type>
	EXPORT constexpr void operator+=(Primitive4<Type>& f1, const Primitive4<Type>& f2)
	{
		f1 = Primitive4<Type>(f1.x + f2.x, f1.y + f2.y, f1.z + f2.z, f1.w + f2.w);
	}
	template <typename Type>
	EXPORT constexpr void operator-=(Primitive2<Type>& f1, const Primitive2<Type>& f2)
	{
		f1 = Primitive2<Type>(f1.x - f2.x, f1.y - f2.y);
	}
	template <typename Type>
	EXPORT constexpr void operator-=(Primitive3<Type>& f1, const Primitive3<Type>& f2)
	{
		f1 = Primitive3<Type>(f1.x - f2.x, f1.y - f2.y, f1.z - f2.z);
	}
	template <typename Type>
	EXPORT constexpr void operator-=(Primitive4<Type>& f1, const Primitive4<Type>& f2)
	{
		f1 = Primitive4<Type>(f1.x - f2.x, f1.y - f2.y, f1.z - f2.z, f1.w - f2.w);
	}
	template <typename Type>
	EXPORT constexpr Primitive2<Type> operator*(const Primitive2<Type>& f1, const Primitive2<Type>& f2)
	{
		return Primitive2<Type>(f1.x * f2.x, f1.y * f2.y);
	}
	template <typename Type>
	EXPORT constexpr Primitive3<Type> operator*(const Primitive3<Type>& f1, const Primitive3<Type>& f2)
	{
		return Primitive3<Type>(f1.x * f2.x, f1.y * f2.y, f1.z * f2.z);
	}
	template <typename Type>
	EXPORT constexpr Primitive4<Type> operator*(const Primitive4<Type>& f1, const Primitive4<Type>& f2)
	{
		return Primitive4<Type>(f1.x * f2.x, f1.y * f2.y, f1.z * f2.z, f1.w * f2.w);
	}
	template <typename Type>
	EXPORT constexpr Primitive2<Type> operator*(const Primitive2<Type>& f1, const double f2)
	{
		return Primitive2<Type>(f1.x * f2, f1.y * f2);
	}
	template <typename Type>
	EXPORT constexpr Primitive3<Type> operator*(const Primitive3<Type>& f1, const double f2)
	{
		return Primitive3<Type>(f1.x * f2, f1.y * f2, f1.z * f2);
	}
	template <typename Type>
	EXPORT constexpr Primitive4<Type> operator*(const Primitive4<Type>& f1, const double f2)
	{
		return Primitive4<Type>(f1.x * f2, f1.y * f2, f1.z * f2, f1.w * f2);
	}

	template <typename Type>
	EXPORT constexpr Primitive2<Type> operator*=(Primitive2<Type>& f1, const Primitive2<Type>& f2)
	{
		f1 = Primitive2<Type>(f1.x * f2.x, f1.y * f2.y);
	}
	template <typename Type>
	EXPORT constexpr Primitive3<Type> operator*=(Primitive3<Type>& f1, const Primitive3<Type>& f2)
	{
		f1 = Primitive3<Type>(f1.x * f2.x, f1.y * f2.y, f1.z * f2.z);
	}
	template <typename Type>
	EXPORT constexpr Primitive4<Type> operator*=(Primitive4<Type>& f1, const Primitive4<Type>& f2)
	{
		f1 = Primitive4<Type>(f1.x * f2.x, f1.y * f2.y, f1.z * f2.z, f1.w * f2.w);
	}
	template <typename Type>
	EXPORT constexpr Primitive2<Type> operator*=(Primitive2<Type>& f1, const double f2)
	{
		f1 = Primitive2<Type>(f1.x * f2, f1.y * f2);
	}
	template <typename Type>
	EXPORT constexpr Primitive3<Type> operator*=(Primitive3<Type>& f1, const double f2)
	{
		f1 = Primitive3<Type>(f1.x * f2, f1.y * f2, f1.z * f2);
	}
	template <typename Type>
	EXPORT constexpr Primitive4<Type> operator*=(Primitive4<Type>& f1, const double f2)
	{
		f1 = Primitive4<Type>(f1.x * f2, f1.y * f2, f1.z * f2, f1.w * f2);
	}
}
#pragma once
#include "mlCommon.hpp"

namespace MathLib
{
	// Primitive types

	class Int4;
	class Int3;
	class Int2;

	class EXPORT Float4
	{
	public:
		union
		{
			struct
			{
				double x;
				double y;
				double z;
				double w;
			};
			struct
			{
				double r;
				double g;
				double b;
				double a;
			};
			double content[4];
		};

		constexpr Float4();
		constexpr Float4(const double x, const double y, const double z, const double w);
		Float4(const Int4& other);
		void operator=(const Int4& other);
		constexpr double operator[](const size_t index);
		constexpr const double operator[](const size_t index) const;
	};

	class EXPORT Float3
	{
	public:
		union
		{
			struct
			{
				double x;
				double y;
				double z;
			};
			struct
			{
				double r;
				double g;
				double b;
			};
			double content[3];
		};

		constexpr Float3();
		constexpr Float3(const double x, const double y, const double z);
		Float3(const Int3& other);
		void operator=(const Int3& other);
		constexpr double operator[](const size_t index);
		constexpr const double operator[](const size_t index) const;
	};

	class EXPORT Float2
	{
	public:
		union
		{
			struct
			{
				double x;
				double y;
			};
			double content[2];
		};

		constexpr Float2();
		constexpr Float2(const double x, const double y);
		Float2(const Int2& other);
		void operator=(const Int2& other);
		constexpr double operator[](const size_t index);
		constexpr const double operator[](const size_t index) const;
	};

	class EXPORT Int4
	{
	public:
		union
		{
			struct
			{
				int x;
				int y;
				int z;
				int w;
			};
			struct
			{
				int r;
				int g;
				int b;
				int a;
			};
			int content[4];
		};

		constexpr Int4()
			: x(0), y(0), z(0), w(0)
		{

		}
		constexpr Int4(const int x, const int y, const int z, const int w)
			: x(x), y(y), z(z), w(w)
		{

		}
		Int4(const Float4& other)
		{
			x = (int) round(other.x);
			y = (int) round(other.y);
			z = (int) round(other.z);
			w = (int) round(other.w);
		}
		void operator=(const Float4& other)
		{
			x = (int) round(other.x);
			y = (int) round(other.y);
			z = (int) round(other.z);
			w = (int) round(other.w);
		}
		constexpr int operator[](const size_t index)
		{
			return content[index];
		}
		constexpr const int operator[](const size_t index) const
		{
			return content[index];
		}
	};

	class EXPORT Int3
	{
	public:
		union
		{
			struct
			{
				int x;
				int y;
				int z;
			};
			struct
			{
				int r;
				int g;
				int b;
			};
			int content[3];
		};

		constexpr Int3()
			: x(0), y(0), z(0)
		{

		}
		constexpr Int3(const int x, const int y, const int z)
			: x(x), y(y), z(z)
		{

		}
		Int3(const Float3& other)
		{
			x = (int) round(other.x);
			y = (int) round(other.y);
			z = (int) round(other.z);
		}
		void operator=(const Float3& other)
		{
			x = (int) round(other.x);
			y = (int) round(other.y);
			z = (int) round(other.z);
		}
		constexpr int operator[](const size_t index)
		{
			return content[index];
		}
		constexpr const int operator[](const size_t index) const
		{
			return content[index];
		}
	};

	class EXPORT Int2
	{
	public:
		union
		{
			struct
			{
				int x;
				int y;
			};
			int content[2];
		};

		constexpr Int2()
			: x(0), y(0)
		{

		}
		constexpr Int2(const int x, const int y)
			: x(x), y(y)
		{

		}
		Int2(const Float2& other)
		{
			x = (int) round(other.x);
			y = (int) round(other.y);
		}
		void operator=(const Float2& other)
		{
			x = (int) round(other.x);
			y = (int) round(other.y);
		}
		constexpr int operator[](const size_t index)
		{
			return content[index];
		}
		constexpr const int operator[](const size_t index) const
		{
			return content[index];
		}
	};

	// Primitives Implementations

	constexpr Float4::Float4()
		: x(0), y(0), z(0), w(0)
	{

	}
	constexpr Float4::Float4(const double x, const double y, const double z, const double w)
		: x(x), y(y), z(z), w(w)
	{

	}
	constexpr double Float4::operator[](const size_t index)
	{
		return content[index];
	}
	constexpr const double Float4::operator[](const size_t index) const
	{
		return content[index];
	}

	constexpr Float3::Float3()
		: x(0), y(0), z(0)
	{

	}
	constexpr Float3::Float3(const double x, const double y, const double z)
		: x(x), y(y), z(z)
	{

	}
	constexpr double Float3::operator[](const size_t index)
	{

		return content[index];

	}
	constexpr const double Float3::operator[](const size_t index) const
	{
		return content[index];
	}

	constexpr Float2::Float2()
		: x(0), y(0)
	{

	}
	constexpr Float2::Float2(const double x, const double y)
		: x(x), y(y)
	{

	}
	constexpr double Float2::operator[](const size_t index)
	{
		return content[index];
	}
	constexpr const double Float2::operator[](const size_t index) const
	{
		return content[index];
	}

	// Operator overloads

	EXPORT constexpr bool operator==(const Float2& f1, const Float2& f2)
	{
		return f1.x == f2.x && f1.y == f2.y;
	}
	EXPORT constexpr bool operator!=(const Float2& f1, const Float2& f2)
	{
		return f1.x != f2.x || f1.y != f2.y;
	}
	EXPORT constexpr bool operator==(const Float3& f1, const Float3& f2)
	{
		return f1.x == f2.x && f1.y == f2.y && f1.z == f2.z;
	}
	EXPORT constexpr bool operator!=(const Float3& f1, const Float3& f2)
	{
		return f1.x != f2.x || f1.y != f2.y || f1.z != f2.z;
	}
	EXPORT constexpr bool operator==(const Float4& f1, const Float4& f2)
	{
		return f1.x == f2.x && f1.y == f2.y && f1.z == f2.z && f1.w == f2.w;
	}
	EXPORT constexpr bool operator!=(const Float4& f1, const Float4& f2)
	{
		return f1.x != f2.x || f1.y != f2.y || f1.z != f2.z || f1.w != f2.w;
	}
	EXPORT constexpr Float2 operator+(const Float2& f1, const Float2& f2)
	{
		return Float2(f1.x + f2.x, f1.y + f2.y);
	}
	EXPORT constexpr Float3 operator+(const Float3& f1, const Float3& f2)
	{
		return Float3(f1.x + f2.x, f1.y + f2.y, f1.z + f2.z);
	}
	EXPORT constexpr Float4 operator+(const Float4& f1, const Float4& f2)
	{
		return Float4(f1.x + f2.x, f1.y + f2.y, f1.z + f2.z, f1.w + f2.w);
	}
	EXPORT constexpr Float2 operator-(const Float2& f1, const Float2& f2)
	{
		return Float2(f1.x - f2.x, f1.y - f2.y);
	}
	EXPORT constexpr Float3 operator-(const Float3& f1, const Float3& f2)
	{
		return Float3(f1.x - f2.x, f1.y - f2.y, f1.z - f2.z);
	}
	EXPORT constexpr Float4 operator-(const Float4& f1, const Float4& f2)
	{
		return Float4(f1.x - f2.x, f1.y - f2.y, f1.z - f2.z, f1.w - f2.w);
	}
	EXPORT constexpr void operator+=(Float2& f1, const Float2& f2)
	{
		f1 = Float2(f1.x + f2.x, f1.y + f2.y);
	}
	EXPORT constexpr void operator+=(Float3& f1, const Float3& f2)
	{
		f1 = Float3(f1.x + f2.x, f1.y + f2.y, f1.z + f2.z);
	}
	EXPORT constexpr void operator+=(Float4& f1, const Float4& f2)
	{
		f1 = Float4(f1.x + f2.x, f1.y + f2.y, f1.z + f2.z, f1.w + f2.w);
	}
	EXPORT constexpr void operator-=(Float2& f1, const Float2& f2)
	{
		f1 = Float2(f1.x - f2.x, f1.y - f2.y);
	}
	EXPORT constexpr void operator-=(Float3& f1, const Float3& f2)
	{
		f1 = Float3(f1.x - f2.x, f1.y - f2.y, f1.z - f2.z);
	}
	EXPORT constexpr void operator-=(Float4& f1, const Float4& f2)
	{
		f1 = Float4(f1.x - f2.x, f1.y - f2.y, f1.z - f2.z, f1.w - f2.w);
	}
	EXPORT constexpr Float2 operator*(const Float2& f1, const Float2& f2)
	{
		return Float2(f1.x * f2.x, f1.y * f2.y);
	}
	EXPORT constexpr Float3 operator*(const Float3& f1, const Float3& f2)
	{
		return Float3(f1.x * f2.x, f1.y * f2.y, f1.z * f2.z);
	}
	EXPORT constexpr Float4 operator*(const Float4& f1, const Float4& f2)
	{
		return Float4(f1.x * f2.x, f1.y * f2.y, f1.z * f2.z, f1.w * f2.w);
	}
	EXPORT constexpr Float2 operator*(const Float2& f1, const double f2)
	{
		return Float2(f1.x * f2, f1.y * f2);
	}
	EXPORT constexpr Float3 operator*(const Float3& f1, const double f2)
	{
		return Float3(f1.x * f2, f1.y * f2, f1.z * f2);
	}
	EXPORT constexpr Float4 operator*(const Float4& f1, const double f2)
	{
		return Float4(f1.x * f2, f1.y * f2, f1.z * f2, f1.w * f2);
	}

	EXPORT constexpr bool operator==(const Int2& i1, const Int2& i2)
	{
		return i1.x == i2.x && i1.y == i2.y;
	}
	EXPORT constexpr bool operator!=(const Int2& i1, const Int2& i2)
	{
		return i1.x != i2.x || i1.y != i2.y;
	}
	EXPORT constexpr bool operator==(const Int3& i1, const Int3& i2)
	{
		return i1.x == i2.x && i1.y == i2.y && i1.z == i2.z;
	}
	EXPORT constexpr bool operator!=(const Int3& i1, const Int3& i2)
	{
		return i1.x != i2.x || i1.y != i2.y || i1.z != i2.z;
	}
	EXPORT constexpr bool operator==(const Int4& i1, const Int4& i2)
	{
		return i1.x == i2.x && i1.y == i2.y && i1.z == i2.z && i1.w == i2.w;
	}
	EXPORT constexpr bool operator!=(const Int4& i1, const Int4& i2)
	{
		return i1.x != i2.x || i1.y != i2.y || i1.z != i2.z || i1.w != i2.w;
	}
	EXPORT constexpr Int2 operator+(const Int2& i1, const Int2& i2)
	{
		return Int2(i1.x + i2.x, i1.y + i2.y);
	}
	EXPORT constexpr Int3 operator+(const Int3& i1, const Int3& i2)
	{
		return Int3(i1.x + i2.x, i1.y + i2.y, i1.z + i2.z);
	}
	EXPORT constexpr Int4 operator+(const Int4& i1, const Int4& i2)
	{
		return Int4(i1.x + i2.x, i1.y + i2.y, i1.z + i2.z, i1.w + i2.w);
	}
	EXPORT constexpr Int2 operator-(const Int2& i1, const Int2& i2)
	{
		return Int2(i1.x - i2.x, i1.y - i2.y);
	}
	EXPORT constexpr Int3 operator-(const Int3& i1, const Int3& i2)
	{
		return Int3(i1.x - i2.x, i1.y - i2.y, i1.z - i2.z);
	}
	EXPORT constexpr Int4 operator-(const Int4& i1, const Int4& i2)
	{
		return Int4(i1.x - i2.x, i1.y - i2.y, i1.z - i2.z, i1.w - i2.w);
	}
	EXPORT constexpr void operator+=(Int2& i1, const Int2& i2)
	{
		i1 = Int2(i1.x + i2.x, i1.y + i2.y);
	}
	EXPORT constexpr void operator+=(Int3& i1, const Int3& i2)
	{
		i1 = Int3(i1.x + i2.x, i1.y + i2.y, i1.z + i2.z);
	}
	EXPORT constexpr void operator+=(Int4& i1, const Int4& i2)
	{
		i1 = Int4(i1.x + i2.x, i1.y + i2.y, i1.z + i2.z, i1.w + i2.w);
	}
	EXPORT constexpr void operator-=(Int2& i1, const Int2& i2)
	{
		i1 = Int2(i1.x - i2.x, i1.y - i2.y);
	}
	EXPORT constexpr void operator-=(Int3& i1, const Int3& i2)
	{
		i1 = Int3(i1.x - i2.x, i1.y - i2.y, i1.z - i2.z);
	}
	EXPORT constexpr void operator-=(Int4& i1, const Int4& i2)
	{
		i1 = Int4(i1.x - i2.x, i1.y - i2.y, i1.z - i2.z, i1.w - i2.w);
	}
	EXPORT constexpr Int2 operator*(const Int2& i1, const Int2& i2)
	{
		return Int2(i1.x * i2.x, i1.y * i2.y);
	}
	EXPORT constexpr Int3 operator*(const Int3& i1, const Int3& i2)
	{
		return Int3(i1.x * i2.x, i1.y * i2.y, i1.z * i2.z);
	}
	EXPORT constexpr Int4 operator*(const Int4& i1, const Int4& i2)
	{
		return Int4(i1.x * i2.x, i1.y * i2.y, i1.z * i2.z, i1.w * i2.w);
	}
	EXPORT constexpr Int2 operator*(const Int2& i1, const int& i2)
	{
		return Int2(i1.x * i2, i1.y * i2);
	}
	EXPORT constexpr Int3 operator*(const Int3& i1, const int& i2)
	{
		return Int3(i1.x * i2, i1.y * i2, i1.z * i2);
	}
	EXPORT constexpr Int4 operator*(const Int4& i1, const int& i2)
	{
		return Int4(i1.x * i2, i1.y * i2, i1.z * i2, i1.w * i2);
	}
}
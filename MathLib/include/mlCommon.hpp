#pragma once

#define MATHLIB

#include <math.h>

#define PI 3.14159265

#ifndef MATHLIB_STATIC
	#ifndef EXPORT
		#ifdef MATHLIB_EXPORTS
			#ifdef _WIN32
				#define EXPORT __declspec(dllexport)
			#elif __unix__
				#define EXPORT __attribute__((visibility("default")))
			#endif // OS check
		#else
			#ifdef _WIN32
				#define EXPORT __declspec(dllimport)
			#endif // __WIN32
		#endif // MATHLIB_EXPORTS
	#endif // EXPORT

#else
	#ifndef EXPORT
		#define EXPORT
	#endif // EXPORT
#endif // MATHLIB_STATIC

#ifndef DEPRECATED
	#if __cplusplus >= 201402L
		#define DEPRECATED(reason) [[deprecated(#reason)]]
	#else
		#define DEPRECATED(reason)
	#endif // __cplusplus version
#endif // DEPRECATED

#ifndef DEPRECATEDCLASS
	#if __cplusplus >= 201402L && defined __WIN32
		#define DEPRECATEDCLASS(reason) __declspec(deprecated(#reason))
	#else
		#define DEPRECATEDCLASS(reason)
	#endif // __cplusplus version
#endif // DEPRECATEDCLASS


namespace MathLib
{
	// Generally useful functions

	template<typename T>
	EXPORT constexpr const T& Max(const T& t1, const T& t2)
	{
		if (t1 > t2)
			return t1;
		return t2;
	}

	template<typename T>
	EXPORT constexpr const T& Min(const T& t1, const T& t2)
	{
		if (t1 < t2)
			return t1;
		return t2;
	}

	template<typename T>
	EXPORT constexpr const T& Clamp(const T& val, const T& min, const T& max)
	{
		if (val > max)
			return max;
		else if (val < min)
			return min;
		return val;
	}

	template<typename T>
	EXPORT constexpr const T& Pow(const T& base, size_t exponent)
	{
		T res = base;
		for (; exponent > 0; exponent--)
			res *= base;
	}
	EXPORT constexpr const size_t PowNeg1(size_t exponent)
	{
		if (exponent % 2 == 0)
			return 1;
		else
			return -1;
	}

	EXPORT constexpr double Deg2Rad(double deg)
	{
		return deg * PI / 180;
	}
	EXPORT constexpr double Rad2Deg(double rad)
	{
		return rad * 180 / PI;
	}
}
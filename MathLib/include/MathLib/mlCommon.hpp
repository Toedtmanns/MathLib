#pragma once

#define MATHLIB

#include <math.h>

#define PI 3.14159265

#ifndef MATHLIB_STATIC
	#ifndef EXPORT
		#ifdef MATHLIB_EXPORTS
			#ifdef _WIN32
				#ifdef _MSC_VER
					#define EXPORT __declspec(dllexport)
				#elif __GNUC__
					#define EXPORT __attribute__((dllexport))
				#endif // OS check
			#else
				#define EXPORT __attribute__((visibility("default")))
			#endif
		#else
			#ifdef _WIN32
				#ifdef _MSC_VER
					#define EXPORT __declspec(dllimport)
				#elif __GNUC__
					#define EXPORT __attribute__((dllimport))
				#endif // __WIN32
			#else
				#define EXPORT
			#endif
		#endif // MATHLIB_EXPORTS
	#endif // EXPORT

#else
	#ifndef EXPORT
		#define EXPORT
	#endif // EXPORT
#endif // MATHLIB_STATIC

#ifndef DEPRECATED
	#if __cplusplus >= 201402L
		#define DEPRECATED(...) [[deprecated(#__VA_ARGS__)]]
	#else
		#define DEPRECATED(...)
	#endif // __cplusplus version
#endif // DEPRECATED

#ifndef DEPRECATEDCLASS
	#if __cplusplus >= 201402L && defined __WIN32
		#define DEPRECATEDCLASS(...) __declspec(deprecated(#__VA_ARGS__))
	#else
		#define DEPRECATEDCLASS(...)
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
	EXPORT constexpr const T& Clamp(const T& val, const T& Min, const T& Max)
	{
		if (val > Max)
			return Max;
		else if (val < Min)
			return Min;
		return val;
	}

	template<typename T>
	EXPORT constexpr const T Pow(const T& base, size_t exponent)
	{
		T res = base;
		for (; exponent > 0; exponent--)
			res *= base;
		return res;
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
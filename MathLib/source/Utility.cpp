#include "../include/MathLib.hpp"
#include <stdio.h>
#include <cmath>
#include <cstring>

namespace MathLib
{
	double Lerp(const double& start, const double& end, const double& t)
	{
		return start + t * (end - start);
	}
	Float2 Lerp(const Float2& start, const Float2& end, const double& t)
	{
		return {
			start.x + t * (end.x - start.x),
			start.y + t * (end.y - start.y)
		};
	}
	Float2 Lerp(const Line2D& line, const double& t)
	{
		return {
			line.p1.x + t * (line.p2.x - line.p1.x),
			line.p1.y + t * (line.p2.y - line.p1.y)
		};
	}

	double MinFromArray(const double* const arr, const unsigned int& length)
	{
		double min = arr[0];
		for (unsigned int i = 1; i < length; i++)
		{
			if (arr[i] < min)
				min = arr[i];
		}
		return min;
	}
	double MaxFromArray(const double* const arr, const unsigned int& length)
	{
		double max = arr[0];
		for (unsigned int i = 1; i < length; i++)
		{
			if (arr[i] > max)
				max = arr[i];
		}
		return max;
	}
	Line2D Vector2Line(Vector2D vector, Float2 pos)
	{
		Line2D line = Line2D();
		line.p1 = pos;
		line.p2.x = vector.x + pos.x;
		line.p2.y = vector.y + pos.y;

		return line;
	}
	int RandInt(int min, int max)
	{
		int range = max + 1 - min;
		return rand() % range + min;
	}
	MatrixF Vec2Mat(const Vector3D& vec, const int& front)
	{
		MatrixF mat = MatrixF(3);
		mat.SetNum(0, 0, vec.x);
		mat.SetNum(1, 1, vec.y);
		mat.SetNum(2, 2, vec.z);
		return mat;
	}
	Vector3D Mat2Vec(const MatrixF& mat, const int& front)
	{
		Vector3D vec = Vector3D();
		vec.x = mat.GetNum(0, 0);
		vec.y = mat.GetNum(1, 1);
		vec.z = mat.GetNum(2, 2);
		return vec;
	}

	bool SetArraySize(void** array, const unsigned int& currLength, const unsigned int& newLength)
	{
		unsigned int newSize = newLength * sizeof(array[0]);

		void* newArray = realloc(*array, newSize);
		if (newArray)
		{
			*array = newArray;
			return true;
		}
		return false;
	}
}
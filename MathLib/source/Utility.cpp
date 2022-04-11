#include "../include/Maths.h"
#include <stdio.h>
#include <cmath>
#include <cstring>

namespace MathLib
{
	double Lerp(const double& start, const double& end, const double& t)
	{
		return start + t * (end - start);
	}
	Primitives::Float2 Lerp(const Primitives::Float2& start, const Primitives::Float2& end, const double& t)
	{
		return {
			start.x + t * (end.x - start.x),
			start.y + t * (end.y - start.y)
		};
	}
	Primitives::Float2 Lerp(const Primitives::Line2D& line, const double& t)
	{
		return {
			line.p1.x + t * (line.p2.x - line.p1.x),
			line.p1.y + t * (line.p2.y - line.p1.y)
		};
	}

	namespace Utility
	{
		double Deg2Rad(double deg)
		{
			return deg * PI / 180;
		}
		double Rad2Deg(double rad)
		{
			return rad * 180 / PI;
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
		Primitives::Line2D Vector2Line(Vectors::Vector2D vector, Primitives::Float2 pos)
		{
			Primitives::Line2D line = Primitives::Line2D();
			line.p1 = pos;
			line.p2.x = vector.direction.x + pos.x;
			line.p2.y = vector.direction.y + pos.y;

			return line;
		}
		Vectors::Vector2D Line2Vector(Primitives::Line2D line)
		{
			Vectors::Vector2D vector = Vectors::Vector2D();
			vector.direction = Primitives::Float2(line.p2.x - line.p1.x, line.p2.y - line.p1.y);
			return vector;
		}
		Vectors::Vector2D Line2Vector(Primitives::Float2 p1, Primitives::Float2 p2)
		{
			Vectors::Vector2D vector = Vectors::Vector2D();
			vector.direction = Primitives::Float2(p2.x - p1.x, p2.y - p1.y);
			return vector;
		}
		int RandInt(int min, int max)
		{
			int range = max + 1 - min;
			return rand() % range + min;
		}
		Matrices::MatrixF Vec2Mat(const Vectors::Vector3D& vec, const int& front)
		{
			Matrices::MatrixF mat = Matrices::MatrixF(3);
			mat.SetNum(0, 0, vec.direction.x);
			mat.SetNum(1, 1, vec.direction.y);
			mat.SetNum(2, 2, vec.direction.z);
			return mat;
		}
		Vectors::Vector3D Mat2Vec(const Matrices::MatrixF& mat, const int& front)
		{
			Vectors::Vector3D vec = Vectors::Vector3D();
			vec.direction.x = mat.GetNum(0, 0);
			vec.direction.y = mat.GetNum(1, 1);
			vec.direction.z = mat.GetNum(2, 2);
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
}
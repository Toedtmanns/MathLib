#include "../include/MathLib/MathLib.hpp"
#include <stdio.h>
#include <cmath>
#include <cstring>

namespace MathLib
{
	float Lerp(const float start, const float end, const float t)
	{
		return start + t * (end - start);
	}
	Float2 Lerp(const Float2& start, const Float2& end, const float t)
	{
		return {
			start.x + t * (end.x - start.x),
			start.y + t * (end.y - start.y)
		};
	}
	Float2 Lerp(const Line2D& line, const float t)
	{
		return {
			line.p1.x + t * (line.p2.x - line.p1.x),
			line.p1.y + t * (line.p2.y - line.p1.y)
		};
	}
	Float3 Lerp(const Float3& start, const Float3& end, const float t)
	{
		return {
			start.x + t * (end.x - start.x),
			start.y + t * (end.y - start.y),
			start.z + t * (end.z - start.z)
		};
	}
	Float3 Lerp(const Line3D& line, const float t)
	{
		return {
			line.p1.x + t * (line.p2.x - line.p1.x),
			line.p1.y + t * (line.p2.y - line.p1.y),
			line.p1.z + t * (line.p2.z - line.p1.z)
		};
	}

	Line2D Vector2Line(Vector2D vector, Float2 pos)
	{
		Line2D line = Line2D();
		line.p1 = pos;
		line.p2.x = vector.x + pos.x;
		line.p2.y = vector.y + pos.y;

		return line;
	}
	int RandInt(int Min, int Max)
	{
		int range = Max + 1 - Min;
		return rand() % range + Min;
	}
}
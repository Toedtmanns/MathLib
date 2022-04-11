#include "../include/Maths.h"
#include <cstring>
#include <stdexcept>

namespace MathLib
{
	namespace Geometry
	{
		GeoCollision::GeoCollision()
			: isColliding(false), m_IntersectArray(new Primitives::Intersect[2]), m_IntersectCount(0), m_IntersectArrLength(2)
		{

		}
		GeoCollision::GeoCollision(GeoCollision&& other) noexcept
			: isColliding(other.isColliding), m_IntersectArray(other.m_IntersectArray), m_IntersectCount(other.m_IntersectCount), m_IntersectArrLength(other.m_IntersectArrLength)
		{
			other.m_IntersectArray = nullptr;
		}
		GeoCollision::GeoCollision(const GeoCollision& other)
			: isColliding(other.isColliding), m_IntersectCount(other.m_IntersectCount), m_IntersectArrLength(other.m_IntersectArrLength), m_IntersectArray(new Primitives::Intersect[other.m_IntersectArrLength])
		{
			memcpy(m_IntersectArray, other.m_IntersectArray, m_IntersectArrLength * sizeof(Primitives::Intersect));
		}
		GeoCollision::GeoCollision(const Primitives::Intersect* intersectArr, const unsigned int& intersectArrLength)
			: m_IntersectCount(intersectArrLength), m_IntersectArrLength(intersectArrLength)
		{
			if (intersectArrLength == 0 || intersectArr == nullptr)
			{
				isColliding = false;
				m_IntersectArray = nullptr;
				m_IntersectCount = 0;
			}
			else
			{
				isColliding = true;
				m_IntersectArray = new Primitives::Intersect[intersectArrLength];
				memcpy(m_IntersectArray, intersectArr, intersectArrLength * sizeof(Primitives::Intersect));
			}
		}
		void GeoCollision::AddIntersect(const Primitives::Intersect& intersect)
		{
			isColliding = true;
			if (m_IntersectCount + 1 >= m_IntersectArrLength)
			{
				if (Utility::SetArraySize((void**) &m_IntersectArray, m_IntersectCount, m_IntersectArrLength + 2))
					m_IntersectArrLength += 2;
			}
			m_IntersectArray[m_IntersectCount] = Primitives::Intersect(intersect);
			m_IntersectCount++;
		}
		unsigned int GeoCollision::GetIntersectCount() const
		{
			return m_IntersectCount;
		}
		void GeoCollision::operator=(GeoCollision&& other) noexcept
		{
			if (m_IntersectArray != nullptr)
				delete[] m_IntersectArray;

			m_IntersectArray = other.m_IntersectArray;
			other.m_IntersectArray = nullptr;
			m_IntersectCount = other.m_IntersectCount;
			isColliding = other.isColliding;
		}
		void GeoCollision::operator=(const GeoCollision& other)
		{
			if (m_IntersectCount != other.m_IntersectCount)
			{
				if (m_IntersectArray != nullptr)
					delete[] m_IntersectArray;
				m_IntersectArray = new Primitives::Intersect[other.m_IntersectCount];
			}
			memcpy(m_IntersectArray, other.m_IntersectArray, other.m_IntersectCount * sizeof(Primitives::Intersect));
			m_IntersectCount = other.m_IntersectCount;
			isColliding = other.isColliding;
		}
		Primitives::Intersect& GeoCollision::operator[](const unsigned int& index)
		{
			if (index < m_IntersectCount)
				return m_IntersectArray[index];
			else
				throw std::out_of_range("Intersect array access out of bounds!");
		}
		GeoCollision::Iterator GeoCollision::begin()
		{
			return GeoCollision::Iterator(m_IntersectArray);
		}
		GeoCollision::Iterator GeoCollision::end()
		{
			return GeoCollision::Iterator(m_IntersectArray + m_IntersectCount);
		}
		GeoCollision::~GeoCollision()
		{
			if (m_IntersectArray != nullptr)
				delete[] m_IntersectArray;
		}

		Polygon2D::Polygon2D()
			: m_PointArr(new Primitives::Float2[3]), m_Corners(3)
		{
			m_PointArr[0] = Primitives::Float2(0, 1);
			m_PointArr[1] = Primitives::Float2(1, -1);
			m_PointArr[2] = Primitives::Float2(-1, -1);
		}
		Polygon2D::Polygon2D(Polygon2D&& other) noexcept
		{
			m_Corners = other.m_Corners;
			m_PointArr = other.m_PointArr;
		}
		Polygon2D::Polygon2D(const Polygon2D& other)
			: m_PointArr(new Primitives::Float2[other.m_Corners]), m_Corners(other.m_Corners)
		{
			memcpy(m_PointArr, other.m_PointArr, m_Corners * sizeof(Primitives::Float2));
		}
		Polygon2D::Polygon2D(const unsigned int& corners)
			: m_Corners(corners), m_PointArr(new Primitives::Float2[corners])
		{

		}
		Polygon2D::Polygon2D(const Primitives::Float2* const pointArr, const unsigned int& corners)
			: m_PointArr(new Primitives::Float2[corners]), m_Corners(corners)
		{
			memcpy(m_PointArr, pointArr, corners * sizeof(double));
		}
		void Polygon2D::Translate(const Primitives::Float2& translation)
		{
			for (unsigned int c = 0; c < m_Corners; c++)
				m_PointArr[c] += translation;
		}
		void Polygon2D::Translate(const double& translateX, const double& translateY)
		{
			for (unsigned int c = 0; c < m_Corners; c++)
				m_PointArr[c] += Primitives::Float2(translateX, translateY);
		}
		void Polygon2D::Rotate(const double& angle)
		{
			Matrices::MatrixF rotMat = Matrices::MatrixF(2, 2);
			rotMat[0][0] = sin(Utility::Deg2Rad(angle + 90));
			rotMat[0][1] = cos(Utility::Deg2Rad(angle + 90));
			rotMat[1][0] = sin(Utility::Deg2Rad(angle));
			rotMat[1][1] = cos(Utility::Deg2Rad(angle));

			Primitives::Float2 center = GetCenter();
			for (unsigned int c = 0; c < m_Corners; c++)
				m_PointArr[c] = (m_PointArr[c] - center) * rotMat + center;
		}
		void Polygon2D::Scale(const double& scale)
		{
			for (Primitives::Float2& corner : *this)
			{
				corner.x = corner.x * scale;
				corner.y = corner.y * scale;
			}
		}
		void Polygon2D::Scale(const double& scaleX, const double& scaleY)
		{
			for (Primitives::Float2& corner : *this)
			{
				corner.x = corner.x * scaleX;
				corner.y = corner.y * scaleY;
			}
		}
		Primitives::Float2 Polygon2D::GetCenter() const
		{
			Primitives::Float2 center = Primitives::Float2();
			for (unsigned int c = 0; c < m_Corners; c++)
				center = center + m_PointArr[c];
			center.x = center.x / m_Corners;
			center.y = center.y / m_Corners;
			return center;
		}
		void Polygon2D::operator=(Polygon2D&& other) noexcept
		{
			delete[] m_PointArr;
			m_Corners = other.m_Corners;
			m_PointArr = other.m_PointArr;
		}
		void Polygon2D::operator=(const Polygon2D& other)
		{
			if (m_Corners == other.m_Corners)
				memcpy(m_PointArr, other.m_PointArr, m_Corners * sizeof(Primitives::Float2));
			else
			{
				delete[] m_PointArr;
				m_Corners = other.m_Corners;
				m_PointArr = new Primitives::Float2[m_Corners];
				memcpy(m_PointArr, other.m_PointArr, m_Corners * sizeof(Primitives::Float2));
			}
		}
		Polygon2D::Iterator Polygon2D::begin()
		{
			return Polygon2D::Iterator(m_PointArr);
		}
		Polygon2D::Iterator Polygon2D::end()
		{
			return Polygon2D::Iterator(m_PointArr + m_Corners);
		}
		Polygon2D::~Polygon2D()
		{
			delete[] m_PointArr;
		}

		Triangle2D::Triangle2D()
			: Polygon2D()
		{

		}
		Triangle2D::Triangle2D(Triangle2D&& other) noexcept
			: Polygon2D(std::move(other))
		{

		}
		Triangle2D::Triangle2D(const Triangle2D& other)
			: Polygon2D(other)
		{

		}
		Triangle2D::Triangle2D(const Primitives::Float2& p1, const Primitives::Float2& p2, const Primitives::Float2& p3)
			: Polygon2D(3)
		{
			m_PointArr[0] = p1;
			m_PointArr[1] = p2;
			m_PointArr[2] = p3;
		}
		Triangle2D::Triangle2D(const Primitives::Float2* const pointArr)
			: Polygon2D(pointArr, 3)
		{
			
		}
		void Triangle2D::operator=(Triangle2D&& other) noexcept
		{
			Polygon2D::operator=(std::move(other));
		}
		void Triangle2D::operator=(const Triangle2D& other)
		{
			Polygon2D::operator=(other);
		}
		bool Triangle2D::CollidesWith(const Triangle2D& other) const
		{
			for (int side = 0; side < 3; side++)
			{
				Vectors::Vector2D viewVec{{
					m_PointArr[side],
					m_PointArr[(side + 1) % 3]
				}};

				double* this1D = ProjectTo1D(*this, viewVec);
				double* other1D = ProjectTo1D(other, viewVec);

				if (!(Utility::MinFromArray(this1D, 3) < Utility::MaxFromArray(other1D, 3) &&
					Utility::MinFromArray(other1D, 3) < Utility::MaxFromArray(this1D, 3)))
				{
					delete[] this1D, other1D;
					return false;
				}
				delete[] this1D, other1D;
			}
			for (int side = 0; side < 3; side++)
			{
				Vectors::Vector2D viewVec{{
						other.m_PointArr[side],
						other.m_PointArr[(side + 1) % 3]
					}};

				double* this1D = ProjectTo1D(*this, viewVec);
				double* other1D = ProjectTo1D(other, viewVec);

				if (!(Utility::MinFromArray(this1D, 3) < Utility::MaxFromArray(other1D, 3) &&
					Utility::MinFromArray(other1D, 3) < Utility::MaxFromArray(this1D, 3)))
				{
					delete[] this1D, other1D;
					return false;
				}
				delete[] this1D, other1D;
			}
			return true;
		}
		GeoCollision Triangle2D::GetCollision(const Triangle2D& other) const
		{
			GeoCollision retCollision = GeoCollision();

			if (this->CollidesWith(other))
			{
				for (unsigned int side = 0; side < m_Corners; side++)
				{
					Primitives::Line2D line = Primitives::Line2D(
						m_PointArr[side],
						m_PointArr[(side + 1) % 3]
					);

					for (unsigned int oSide = 0; oSide < other.m_Corners; oSide++)
					{
						Primitives::Line2D oLine = Primitives::Line2D(
							other.m_PointArr[oSide],
							other.m_PointArr[(oSide + 1) % 3]
						);

						Primitives::Intersect intersect = Primitives::GetIntersect(line, oLine);
						if (intersect.isIntersecting)
							retCollision.AddIntersect(intersect);
					}
				}
			}

			return retCollision;
		}

		Rectangle2D::Rectangle2D()
			: Polygon2D(4)
		{
			m_PointArr[0] = Primitives::Float2(1, 1);
			m_PointArr[1] = Primitives::Float2(1, -1);
			m_PointArr[2] = Primitives::Float2(-1, -1);
			m_PointArr[3] = Primitives::Float2(-1, 1);
		}
		Rectangle2D::Rectangle2D(Rectangle2D&& other) noexcept
			: Polygon2D(std::move(other))
		{

		}
		Rectangle2D::Rectangle2D(const Rectangle2D& other)
			: Polygon2D(other)
		{

		}
		Rectangle2D::Rectangle2D(const Primitives::Float2& p1, const Primitives::Float2& p2, const Primitives::Float2& p3, const Primitives::Float2& p4)
			: Polygon2D(4)
		{
			m_PointArr[0] = p1;
			m_PointArr[1] = p2;
			m_PointArr[2] = p3;
			m_PointArr[3] = p4;
		}
		Rectangle2D::Rectangle2D(const Primitives::Float2* const pointArr)
			: Polygon2D(pointArr, 4)
		{

		}
		Rectangle2D::Rectangle2D(const Primitives::Float2 position, const double& rotation, const Primitives::Float2 scale)
			: Polygon2D(4)
		{
			m_PointArr[0] = Primitives::Float2(1, 1);
			m_PointArr[1] = Primitives::Float2(1, -1);
			m_PointArr[2] = Primitives::Float2(-1, -1);
			m_PointArr[3] = Primitives::Float2(-1, 1);

			Scale(scale.x, scale.y);
			Rotate(rotation);
			Translate(position);
		}
		void Rectangle2D::operator=(Rectangle2D&& other) noexcept
		{
			Polygon2D::operator=(std::move(other));
		}
		void Rectangle2D::operator=(const Rectangle2D& other)
		{
			Polygon2D::operator=(other);
		}
		bool Rectangle2D::CollidesWith(const Rectangle2D& other) const
		{
			for (int side = 0; side < 2; side++)
			{
				Vectors::Vector2D viewVec{{
						m_PointArr[side],
						m_PointArr[(side + 1) % 2]
					}};

				double* this1D = ProjectTo1D(*this, viewVec);
				double* other1D = ProjectTo1D(other, viewVec);

				if (!(Utility::MinFromArray(this1D, 4) < Utility::MaxFromArray(other1D, 4) &&
					Utility::MinFromArray(other1D, 4) < Utility::MaxFromArray(this1D, 4)))
				{
					delete[] this1D, other1D;
					return false;
				}
				delete[] this1D, other1D;
			}
			for (int side = 0; side < 2; side++)
			{
				Vectors::Vector2D viewVec{{
						other.m_PointArr[side],
						other.m_PointArr[(side + 1) % 2]
					}};

				double* this1D = ProjectTo1D(*this, viewVec);
				double* other1D = ProjectTo1D(other, viewVec);

				if (!(Utility::MinFromArray(this1D, 4) < Utility::MaxFromArray(other1D, 4) &&
					Utility::MinFromArray(other1D, 4) < Utility::MaxFromArray(this1D, 4)))
				{
					delete[] this1D, other1D;
					return false;
				}
				delete[] this1D, other1D;
			}
			return true;
		}
		GeoCollision Rectangle2D::GetCollision(const Rectangle2D& other) const
		{
			GeoCollision retCollision = GeoCollision();

			if (this->CollidesWith(other))
			{
				for (unsigned int side = 0; side < m_Corners; side++)
				{
					Primitives::Line2D line = Primitives::Line2D(
						m_PointArr[side],
						m_PointArr[(side + 1) % 4]
					);

					for (unsigned int oSide = 0; oSide < other.m_Corners; oSide++)
					{
						Primitives::Line2D oLine = Primitives::Line2D(
							other.m_PointArr[oSide],
							other.m_PointArr[(oSide + 1) % 4]
						);

						Primitives::Intersect intersect = Primitives::GetIntersect(line, oLine);
						if (intersect.isIntersecting)
							retCollision.AddIntersect(intersect);
					}
				}
			}

			return retCollision;
		}

		bool Contains(const Rectangle2D& rect, const double& rotation, const Primitives::Float2& point)
		{
			Rectangle2D calcRect = rect;
			calcRect.Rotate(-rotation);
			Vectors::Vector2D pointVec = {Primitives::Line2D(rect.GetCenter(), point)};
			pointVec.Rotate(-rotation);
			Primitives::Float2 transPoint = pointVec.TransformPoint(calcRect.GetCenter());

			double xArr[2] = {calcRect.m_PointArr[0].x, calcRect.m_PointArr[2].x};
			double yArr[2] = {calcRect.m_PointArr[0].y, calcRect.m_PointArr[2].y};

			return transPoint.x < Utility::MaxFromArray(xArr, 2) && transPoint.x > Utility::MinFromArray(xArr, 2) &&
				transPoint.y < Utility::MaxFromArray(yArr, 2) && transPoint.y > Utility::MinFromArray(yArr, 2);
		}
		double* ProjectTo1D(const Polygon2D& polygon, const Vectors::Vector2D& viewDir)
		{
			double* retArr = new double[polygon.m_Corners];

			for (unsigned int p = 0; p < polygon.m_Corners; p++)
			{
				Vectors::Vector2D transVec = Vectors::Vector2D(polygon.m_PointArr[p]);
				transVec.Rotate(-(viewDir.GetAngle()));
				retArr[p] = -transVec.direction.y;
			}

			return retArr;
		}
		double* ProjectTo1D(const Primitives::Float2* pointArr, const unsigned int& pointCount, const Vectors::Vector2D& viewDir)
		{
			double* retArr = new double[pointCount];

			for (unsigned int p = 0; p < pointCount; p++)
			{
				Vectors::Vector2D transVec = Vectors::Vector2D(pointArr[p]);
				transVec.Rotate(-(viewDir.GetAngle()));
				retArr[p] = -transVec.direction.y;
			}

			return retArr;
		}
	}
}
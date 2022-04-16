#include "../include/Maths.hpp"
#include <cstring>
#include <stdexcept>

namespace MathLib
{
	GeoCollision::GeoCollision()
		: isColliding(false), m_IntersectArray(new Intersect[2]), m_IntersectCount(0), m_IntersectArrLength(2)
	{

	}
	GeoCollision::GeoCollision(GeoCollision&& other) noexcept
		: isColliding(other.isColliding), m_IntersectArray(other.m_IntersectArray), m_IntersectCount(other.m_IntersectCount), m_IntersectArrLength(other.m_IntersectArrLength)
	{
		other.m_IntersectArray = nullptr;
	}
	GeoCollision::GeoCollision(const GeoCollision& other)
		: isColliding(other.isColliding), m_IntersectCount(other.m_IntersectCount), m_IntersectArrLength(other.m_IntersectArrLength), m_IntersectArray(new Intersect[other.m_IntersectArrLength])
	{
		memcpy(m_IntersectArray, other.m_IntersectArray, m_IntersectArrLength * sizeof(Intersect));
	}
	GeoCollision::GeoCollision(const Intersect* intersectArr, const unsigned int& intersectArrLength)
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
			m_IntersectArray = new Intersect[intersectArrLength];
			memcpy(m_IntersectArray, intersectArr, intersectArrLength * sizeof(Intersect));
		}
	}
	void GeoCollision::AddIntersect(const Intersect& intersect)
	{
		isColliding = true;
		if (m_IntersectCount + 1 >= m_IntersectArrLength)
		{
			if (SetArraySize((void**) &m_IntersectArray, m_IntersectCount, m_IntersectArrLength + 2))
				m_IntersectArrLength += 2;
		}
		m_IntersectArray[m_IntersectCount] = Intersect(intersect);
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
			m_IntersectArray = new Intersect[other.m_IntersectCount];
		}
		memcpy(m_IntersectArray, other.m_IntersectArray, other.m_IntersectCount * sizeof(Intersect));
		m_IntersectCount = other.m_IntersectCount;
		isColliding = other.isColliding;
	}
	Intersect& GeoCollision::operator[](const unsigned int& index)
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
		: m_CornerArr(new Float2[3]), m_NumCorners(3)
	{
		m_CornerArr[0] = Float2(0, 1);
		m_CornerArr[1] = Float2(1, -1);
		m_CornerArr[2] = Float2(-1, -1);
	}
	Polygon2D::Polygon2D(Polygon2D&& other) noexcept
	{
		m_NumCorners = other.m_NumCorners;
		m_CornerArr = other.m_CornerArr;
	}
	Polygon2D::Polygon2D(const Polygon2D& other)
		: m_CornerArr(new Float2[other.m_NumCorners]), m_NumCorners(other.m_NumCorners)
	{
		memcpy(m_CornerArr, other.m_CornerArr, m_NumCorners * sizeof(Float2));
	}
	Polygon2D::Polygon2D(const unsigned int& corners)
		: m_NumCorners(corners), m_CornerArr(new Float2[corners])
	{

	}
	Polygon2D::Polygon2D(const Float2* const pointArr, const unsigned int& corners)
		: m_CornerArr(new Float2[corners]), m_NumCorners(corners)
	{
		memcpy(m_CornerArr, pointArr, corners * sizeof(Float2));
	}
	void Polygon2D::Translate(const Float2& translation)
	{
		for (unsigned int c = 0; c < m_NumCorners; c++)
			m_CornerArr[c] += translation;
	}
	void Polygon2D::Translate(const double& translateX, const double& translateY)
	{
		for (unsigned int c = 0; c < m_NumCorners; c++)
			m_CornerArr[c] += Float2(translateX, translateY);
	}
	void Polygon2D::Rotate(const double& angle)
	{
		MatrixF rotMat = MatrixF(2, 2);
		rotMat[0][0] = sin(Deg2Rad(angle + 90));
		rotMat[0][1] = cos(Deg2Rad(angle + 90));
		rotMat[1][0] = sin(Deg2Rad(angle));
		rotMat[1][1] = cos(Deg2Rad(angle));

		Float2 center = GetCenter();
		for (unsigned int c = 0; c < m_NumCorners; c++)
			m_CornerArr[c] = (m_CornerArr[c] - center) * rotMat + center;
	}
	void Polygon2D::Scale(const double& scale)
	{
		for (Float2& corner : *this)
		{
			corner.x = corner.x * scale;
			corner.y = corner.y * scale;
		}
	}
	void Polygon2D::Scale(const double& scaleX, const double& scaleY)
	{
		for (Float2& corner : *this)
		{
			corner.x = corner.x * scaleX;
			corner.y = corner.y * scaleY;
		}
	}
	Float2 Polygon2D::GetCenter() const
	{
		Float2 center = Float2();
		for (unsigned int c = 0; c < m_NumCorners; c++)
			center = center + m_CornerArr[c];
		center.x = center.x / m_NumCorners;
		center.y = center.y / m_NumCorners;
		return center;
	}
	Float2 Polygon2D::SupportFunction(const Vector2D& direction) const
	{
		Float2 center = GetCenter();
		double maxDot;
		Float2 matching = m_CornerArr[0];

		maxDot = Vector2D(center, m_CornerArr[0]) * direction;
		for (unsigned int c = 1; c < m_NumCorners; c++)
		{
			double tempVal = Vector2D(center, m_CornerArr[c]) * direction;
			if (tempVal > maxDot)
			{
				maxDot = tempVal;
				matching = m_CornerArr[c];
			}
		}
		return matching;
	}
	bool Polygon2D::CollidesWith(const Polygon2D& other) const
	{
		Vector2D dirThis{1, 1};
		Vector2D originCheck;
		Float2 checkTri[3];

		Vector2D vecAB, vecAC;

		checkTri[0] = SupportFunction(dirThis) - other.SupportFunction(-dirThis);

		dirThis = Vector2D(checkTri[0], {0, 0});

		checkTri[1] = SupportFunction(dirThis) - other.SupportFunction(-dirThis);

		originCheck = Vector2D(checkTri[0], {0, 0});
		if (originCheck * Vector2D({0, 0}, checkTri[1]) < 0)
			return false;

		dirThis = VectorTripleProduct(Vector2D(checkTri[0], checkTri[1]), originCheck);

		int a = 2, b = 1, c = 0;
		while (true)
		{
			checkTri[a] = SupportFunction(dirThis) - other.SupportFunction(-dirThis);

			originCheck = Vector2D(checkTri[c], {0, 0});
			if (dirThis * Vector2D({0, 0}, checkTri[a]) <= 0)
				return false;

			vecAB = Vector2D(checkTri[a], checkTri[b]);
			vecAC = Vector2D(checkTri[a], checkTri[c]);
			originCheck = Vector2D(checkTri[a], {0, 0});

			a = (a + 1) % 3;
			b = (b + 1) % 3;
			c = (c + 1) % 3;

			if ((dirThis = -VectorTripleProduct(vecAB, vecAC)) * originCheck > 0)
				continue;
			else if ((dirThis = -VectorTripleProduct(vecAC, vecAB)) * originCheck > 0)
				continue;
			else
				return true;
		}
	}
	bool Polygon2D::CollidesWith(const Circle2D& other) const
	{
		Vector2D dirThis{1, 1};
		Vector2D originCheck;
		Float2 checkTri[3];

		Vector2D vecAB, vecAC;

		checkTri[0] = SupportFunction(dirThis) - other.SupportFunction(-dirThis);

		dirThis = Vector2D(checkTri[0], {0, 0});

		checkTri[1] = SupportFunction(dirThis) - other.SupportFunction(-dirThis);

		originCheck = Vector2D(checkTri[0], {0, 0});
		if (originCheck * Vector2D({0, 0}, checkTri[1]) < 0)
			return false;

		dirThis = VectorTripleProduct(Vector2D(checkTri[0], checkTri[1]), originCheck);

		int a = 2, b = 1, c = 0;
		while (true)
		{
			checkTri[a] = SupportFunction(dirThis) - other.SupportFunction(-dirThis);

			originCheck = Vector2D(checkTri[c], {0, 0});
			if (dirThis * Vector2D({0, 0}, checkTri[a]) <= 0)
				return false;

			vecAB = Vector2D(checkTri[a], checkTri[b]);
			vecAC = Vector2D(checkTri[a], checkTri[c]);
			originCheck = Vector2D(checkTri[a], {0, 0});

			a = (a + 1) % 3;
			b = (b + 1) % 3;
			c = (c + 1) % 3;

			if ((dirThis = -VectorTripleProduct(vecAB, vecAC)) * originCheck > 0)
				continue;
			else if ((dirThis = -VectorTripleProduct(vecAC, vecAB)) * originCheck > 0)
				continue;
			else
				return true;
		}
	}
	bool Polygon2D::Contains(const Float2& point) const
	{
		Vector2D pointVec = Vector2D(m_CornerArr[0], point);
		Vector2D checkVec = Vector2D(m_CornerArr[0], m_CornerArr[1]);
		checkVec.Rot90R();
		bool neg = false;
		if (pointVec * checkVec < 0)
			neg = true;

		for (unsigned int c = 1; c < m_NumCorners; c++)
		{
			pointVec = Vector2D(m_CornerArr[c], point);
			checkVec = Vector2D(m_CornerArr[c], m_CornerArr[(c + 1) % m_NumCorners]);
			checkVec.Rot90R();

			if (!neg && pointVec * checkVec < 0 || neg && pointVec * checkVec >= 0)
				return false;
		}
		return true;
	}
	void Polygon2D::operator=(Polygon2D&& other) noexcept
	{
		delete[] m_CornerArr;
		m_NumCorners = other.m_NumCorners;
		m_CornerArr = other.m_CornerArr;
	}
	void Polygon2D::operator=(const Polygon2D& other)
	{
		if (m_NumCorners == other.m_NumCorners)
			memcpy(m_CornerArr, other.m_CornerArr, m_NumCorners * sizeof(Float2));
		else
		{
			delete[] m_CornerArr;
			m_NumCorners = other.m_NumCorners;
			m_CornerArr = new Float2[m_NumCorners];
			memcpy(m_CornerArr, other.m_CornerArr, m_NumCorners * sizeof(Float2));
		}
	}
	Polygon2D::Iterator Polygon2D::begin()
	{
		return Polygon2D::Iterator(m_CornerArr);
	}
	Polygon2D::Iterator Polygon2D::end()
	{
		return Polygon2D::Iterator(m_CornerArr + m_NumCorners);
	}
	Polygon2D::~Polygon2D()
	{
		delete[] m_CornerArr;
	}

	// Circle2D class definitions

	Circle2D::Circle2D()
		: position(), radius(0.5)
	{

	}
	Circle2D::Circle2D(const Float2& position, const double& radius)
		: position(position), radius(radius)
	{

	}
	void Circle2D::Translate(const Float2& translation)
	{
		position += translation;
	}
	void Circle2D::Translate(const double& translateX, const double& translateY)
	{
		position.x += translateX;
		position.y += translateY;
	}
	Float2 Circle2D::SupportFunction(const Vector2D& direction) const
	{
		Vector2D calcVec = direction;
		calcVec.SetScale(radius);
		return (Float2) (calcVec) + position;
	}
	bool Circle2D::CollidesWith(const Polygon2D& other) const
	{
		return other.CollidesWith(*this);
	}
	bool Circle2D::CollidesWith(const Circle2D& other) const
	{
		return GetDistance(position, other.position) < radius + other.radius;
	}

	// Triangle2D class definitions

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
	Triangle2D::Triangle2D(const Float2& p1, const Float2& p2, const Float2& p3)
		: Polygon2D(3)
	{
		m_CornerArr[0] = p1;
		m_CornerArr[1] = p2;
		m_CornerArr[2] = p3;
	}
	Triangle2D::Triangle2D(const Float2* const pointArr)
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
			Vector2D viewVec{{
					m_CornerArr[side],
					m_CornerArr[(side + 1) % 3]
				}};

			double* this1D = ProjectTo1D(*this, viewVec);
			double* other1D = ProjectTo1D(other, viewVec);

			if (!(MinFromArray(this1D, 3) < MaxFromArray(other1D, 3) &&
				MinFromArray(other1D, 3) < MaxFromArray(this1D, 3)))
			{
				delete[] this1D, other1D;
				return false;
			}
			delete[] this1D, other1D;
		}
		for (int side = 0; side < 3; side++)
		{
			Vector2D viewVec{{
					other.m_CornerArr[side],
					other.m_CornerArr[(side + 1) % 3]
				}};

			double* this1D = ProjectTo1D(*this, viewVec);
			double* other1D = ProjectTo1D(other, viewVec);

			if (!(MinFromArray(this1D, 3) < MaxFromArray(other1D, 3) &&
				MinFromArray(other1D, 3) < MaxFromArray(this1D, 3)))
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
			for (unsigned int side = 0; side < m_NumCorners; side++)
			{
				Line2D line = Line2D(
					m_CornerArr[side],
					m_CornerArr[(side + 1) % 3]
				);

				for (unsigned int oSide = 0; oSide < other.m_NumCorners; oSide++)
				{
					Line2D oLine = Line2D(
						other.m_CornerArr[oSide],
						other.m_CornerArr[(oSide + 1) % 3]
					);

					Intersect intersect = GetIntersect(line, oLine);
					if (intersect.isIntersecting)
						retCollision.AddIntersect(intersect);
				}
			}
		}

		return retCollision;
	}

	// Rectangle2D class definitions

	Rectangle2D::Rectangle2D()
		: Polygon2D(4)
	{
		m_CornerArr[0] = Float2(0.5, 0.5);
		m_CornerArr[1] = Float2(0.5, -0.5);
		m_CornerArr[2] = Float2(-0.5, -0.5);
		m_CornerArr[3] = Float2(-0.5, 0.5);
	}
	Rectangle2D::Rectangle2D(Rectangle2D&& other) noexcept
		: Polygon2D(std::move(other))
	{

	}
	Rectangle2D::Rectangle2D(const Rectangle2D& other)
		: Polygon2D(other)
	{

	}
	Rectangle2D::Rectangle2D(const Float2& p1, const Float2& p2, const Float2& p3, const Float2& p4)
		: Polygon2D(4)
	{
		m_CornerArr[0] = p1;
		m_CornerArr[1] = p2;
		m_CornerArr[2] = p3;
		m_CornerArr[3] = p4;
	}
	Rectangle2D::Rectangle2D(const Float2* const pointArr)
		: Polygon2D(pointArr, 4)
	{

	}
	Rectangle2D::Rectangle2D(const Float2 position, const double& rotation, const Float2 scale)
		: Polygon2D(4)
	{
		m_CornerArr[0] = Float2(0.5, 0.5);
		m_CornerArr[1] = Float2(0.5, -0.5);
		m_CornerArr[2] = Float2(-0.5, -0.5);
		m_CornerArr[3] = Float2(-0.5, 0.5);

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
			Vector2D viewVec{{
					m_CornerArr[side],
					m_CornerArr[(side + 1) % 2]
				}};

			double* this1D = ProjectTo1D(*this, viewVec);
			double* other1D = ProjectTo1D(other, viewVec);

			if (!(MinFromArray(this1D, 4) < MaxFromArray(other1D, 4) &&
				MinFromArray(other1D, 4) < MaxFromArray(this1D, 4)))
			{
				delete[] this1D, other1D;
				return false;
			}
			delete[] this1D, other1D;
		}
		for (int side = 0; side < 2; side++)
		{
			Vector2D viewVec{{
					other.m_CornerArr[side],
					other.m_CornerArr[(side + 1) % 2]
				}};

			double* this1D = ProjectTo1D(*this, viewVec);
			double* other1D = ProjectTo1D(other, viewVec);

			if (!(MinFromArray(this1D, 4) < MaxFromArray(other1D, 4) &&
				MinFromArray(other1D, 4) < MaxFromArray(this1D, 4)))
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
			for (unsigned int side = 0; side < m_NumCorners; side++)
			{
				Line2D line = Line2D(
					m_CornerArr[side],
					m_CornerArr[(side + 1) % 4]
				);

				for (unsigned int oSide = 0; oSide < other.m_NumCorners; oSide++)
				{
					Line2D oLine = Line2D(
						other.m_CornerArr[oSide],
						other.m_CornerArr[(oSide + 1) % 4]
					);

					Intersect intersect = GetIntersect(line, oLine);
					if (intersect.isIntersecting)
						retCollision.AddIntersect(intersect);
				}
			}
		}

		return retCollision;
	}

	bool Contains(const Rectangle2D& rect, const double& rotation, const Float2& point)
	{
		Rectangle2D calcRect = rect;
		calcRect.Rotate(-rotation);
		Vector2D pointVec = {Line2D(rect.GetCenter(), point)};
		pointVec.Rotate(-rotation);
		Float2 transPoint = pointVec.TransformPoint(calcRect.GetCenter());

		double xArr[2] = {calcRect.m_CornerArr[0].x, calcRect.m_CornerArr[2].x};
		double yArr[2] = {calcRect.m_CornerArr[0].y, calcRect.m_CornerArr[2].y};

		return transPoint.x < MaxFromArray(xArr, 2) && transPoint.x > MinFromArray(xArr, 2) &&
			transPoint.y < MaxFromArray(yArr, 2) && transPoint.y > MinFromArray(yArr, 2);
	}
	double* ProjectTo1D(const Polygon2D& polygon, const Vector2D& viewDir)
	{
		double* retArr = new double[polygon.m_NumCorners];

		for (unsigned int p = 0; p < polygon.m_NumCorners; p++)
		{
			Vector2D transVec = Vector2D(polygon.m_CornerArr[p]);
			transVec.Rotate(-(viewDir.GetAngle()));
			retArr[p] = -transVec.y;
		}

		return retArr;
	}
	double* ProjectTo1D(const Float2* pointArr, const unsigned int& pointCount, const Vector2D& viewDir)
	{
		double* retArr = new double[pointCount];

		for (unsigned int p = 0; p < pointCount; p++)
		{
			Vector2D transVec = Vector2D(pointArr[p]);
			transVec.Rotate(-(viewDir.GetAngle()));
			retArr[p] = -transVec.y;
		}

		return retArr;
	}
}
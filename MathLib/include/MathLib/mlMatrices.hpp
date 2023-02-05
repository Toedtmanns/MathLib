#pragma once
#include "mlPrimitives.hpp"
#include <string.h>
#include <limits.h>

namespace MathLib
{
	// Matrix maths

	template <size_t columns, size_t rows = columns>
	class Matrix
	{
	protected:
		float m_Matrix[columns * rows]{0};

	public:
		constexpr Matrix(Matrix<columns, rows>&& other) noexcept
		{
			memcpy(m_Matrix, other.m_Matrix, sizeof(m_Matrix));
		}
		constexpr Matrix(const Matrix<columns, rows>& other)
		{
			memcpy(m_Matrix, other.m_Matrix, sizeof(m_Matrix));
		}
		template <size_t oColumns, size_t oRows>
		constexpr Matrix(const Matrix<oColumns, oRows>& other)
		{
			for (size_t i = 0, j = 0; i < columns * rows && j < oColumns * oRows; i += rows, j += oRows)
			{
				memcpy(m_Matrix + i, other.GetMatrix() + j, Min(rows, oRows) * sizeof(float));
			}
		}
		constexpr Matrix()
		{

		}
		explicit constexpr Matrix(const float defaultVal)
		{
			if (columns != rows)
				return;

			for (size_t i = 0; i < columns * rows; i += columns + 1)
			{
				m_Matrix[i] = defaultVal;
			}
		}

		constexpr Matrix<columns, rows>& SetVal(const size_t column, const size_t row, const float value)
		{
			m_Matrix[column * rows + row] = value;
			return *this;
		}
		constexpr Matrix<columns, rows>& SetVal(const size_t index, const float value)
		{
			m_Matrix[index] = value;
			return *this;
		}
		constexpr Matrix<columns, rows>& SetMatrix(const float matValues[columns * rows])
		{
			memcpy(m_Matrix, matValues, sizeof(m_Matrix));
			return *this;
		}
		constexpr Matrix<columns, rows>& SetRow(const size_t row, const float* rowData)
		{
			for (size_t c = 0; c < columns; c++)
			{
				m_Matrix[c * rows + row] = rowData[c];
			}
			return *this;
		}
		constexpr Matrix<columns, rows>& SetRow(const size_t row, const Matrix<columns, 1>& rowData)
		{
			for (size_t c = 0; c < columns; c++)
			{
				m_Matrix[c * rows + row] = rowData.GetVal(c);
			}
			return *this;
		}
		constexpr Matrix<columns, rows>& MultRow(const size_t row, const size_t value)
		{
			for (size_t c = 0; c < columns; c++)
			{
				m_Matrix[c * rows + row] *= value;
			}
			return *this;
		}

		// Getter methods

		constexpr Matrix<columns, 1> GetRow(const size_t row) const
		{
			Matrix<columns, 1> rowMat;
			for (size_t c = 0; c < columns; c++)
			{
				rowMat.SetVal(c, m_Matrix[c * rows + row]);
			}
			return rowMat;
		}
		constexpr Matrix<1, rows> GetColumn(const size_t col) const
		{
			Matrix<1, columns> colMat;
			for (size_t r = 0; r < rows; r++)
			{
				colMat.SetVal(r, m_Matrix[col * rows + r]);
			}
			return colMat;
		}
		constexpr float GetVal(const size_t column, const size_t row) const
		{
			return m_Matrix[column * rows + row];
		}
		constexpr float GetVal(const size_t index) const
		{
			return m_Matrix[index];
		}
		constexpr float* GetMatrix()
		{
			return m_Matrix;
		}
		constexpr const float* GetMatrix() const
		{
			return m_Matrix;
		}
		constexpr const size_t GetColumnCount() const
		{
			return columns;
		}
		constexpr const size_t GetRowCount() const
		{
			return rows;
		}

		// Operator overloads

		void operator=(const Matrix<columns, rows>& other)
		{
			memcpy(m_Matrix, other.m_Matrix, sizeof(m_Matrix));
		}
		void operator=(Matrix<columns, rows>&& other) noexcept
		{
			memcpy(m_Matrix, other.m_Matrix, sizeof(m_Matrix));
		}
		constexpr const float* operator[](const size_t index) const
		{
			return m_Matrix + index * rows;
		}
		constexpr float* operator[](const size_t index)
		{
			return m_Matrix + index * rows;
		}

		constexpr Matrix<columns, rows> operator*(const float value)
		{
			Matrix<columns, rows> retMat{0};
			for (size_t i = 0; i < columns * rows; i++)
			{
				retMat.SetVal(i, m_Matrix[i] * value);
			}
			return retMat;
		}
		constexpr void operator*=(const float value)
		{
			for (size_t i = 0; i < columns * rows; i++)
			{
				m_Matrix[i] = m_Matrix[i] * value;
			}
		}
		template<size_t oColumns>
		constexpr Matrix<oColumns, rows> operator*(const Matrix<oColumns, columns>& other) const
		{
			Matrix<oColumns, rows> retMat;
			for (size_t c = 0; c < oColumns; c++)
			{
				for (size_t r = 0; r < rows; r++)
				{
					float num = 0;
					for (size_t calcCol = 0; calcCol < columns; calcCol++)
					{
						num = num + m_Matrix[calcCol * rows + r] * other.GetMatrix()[c * columns + calcCol];
					}
					retMat.SetVal(c, r, num);
				}
			}
			return retMat;
		}
		template<size_t oColumns>
		constexpr void operator*=(const Matrix<oColumns, columns>& other)
		{
			for (size_t c = 0; c < oColumns; c++)
			{
				for (size_t r = 0; r < rows; r++)
				{
					float num = 0;
					for (size_t calcCol = 0; calcCol < columns; calcCol++)
					{
						num = num + m_Matrix[calcCol * rows + r] * other.GetMatrix()[c * columns + calcCol];
					}
					m_Matrix[c * rows + r] = num;
				}
			}
		}
		constexpr Matrix<columns, rows> operator+(const float value) const
		{
			Matrix<columns, rows> retMat;
			for (size_t i = 0; i < columns * rows; i++)
				retMat.SetVal(i, m_Matrix[i] + value);
			return retMat;
		}
		constexpr void operator+=(const float value)
		{
			for (size_t i = 0; i < columns * rows; i++)
				m_Matrix[i] += value;
		}
		constexpr Matrix<columns, rows> operator+(const Matrix<columns, rows>& other) const
		{
			Matrix<columns, rows> retMat;
			for (size_t i = 0; i < columns * rows; i++)
				retMat.SetVal(i, m_Matrix[i] + other.GetVal(i));
			return retMat;
		}
		constexpr void operator+=(const Matrix<columns, rows>& other)
		{
			for (size_t i = 0; i < columns * rows; i++)
				m_Matrix[i] += other.GetVal(i);
		}
		constexpr Matrix<columns, rows> operator-(const float value) const
		{
			Matrix<columns, rows> retMat;
			for (size_t i = 0; i < columns * rows; i++)
				retMat.SetVal(i, m_Matrix[i] - value);
			return retMat;
		}
		constexpr void operator-=(const float value)
		{
			for (size_t i = 0; i < columns * rows; i++)
				m_Matrix[i] -= value;
		}
		constexpr Matrix<columns, rows> operator-(const Matrix<columns, rows>& other) const
		{
			Matrix<columns, rows> retMat;
			for (size_t i = 0; i < columns * rows; i++)
				retMat.SetVal(i, m_Matrix[i] - other.GetVal(i));
			return retMat;
		}
		constexpr void operator-=(const Matrix<columns, rows>& other)
		{
			for (size_t i = 0; i < columns * rows; i++)
				m_Matrix[i] -= other.GetVal(i);
		}
	};

	typedef Matrix<4> Mat4;
	typedef Matrix<3> Mat3;
	typedef Matrix<2> Mat2;
	typedef Matrix<1, 2> Mat1x2;
	typedef Matrix<1, 3> Mat1x3;
	typedef Matrix<1, 4> Mat1x4;
	typedef Matrix<2, 1> Mat2x1;
	typedef Matrix<2, 3> Mat2x3;
	typedef Matrix<2, 4> Mat2x4;
	typedef Matrix<3, 1> Mat3x1;
	typedef Matrix<3, 2> Mat3x2;
	typedef Matrix<3, 4> Mat3x4;
	typedef Matrix<4, 1> Mat4x1;
	typedef Matrix<4, 2> Mat4x2;
	typedef Matrix<4, 3> Mat4x3;

	constexpr Float2 operator*(const Float2& point, const Mat2& matrix)
	{
		return Float2{
			point.x * matrix.GetVal(0) + point.y * matrix.GetVal(1),
			point.x * matrix.GetVal(2) + point.y * matrix.GetVal(3)
		};
	}
	constexpr void operator*=(Float2& point, const Mat2& matrix)
	{
		point.x = point.x * matrix.GetVal(0) + point.y * matrix.GetVal(1);
		point.y = point.x * matrix.GetVal(2) + point.y * matrix.GetVal(3);
	}
	constexpr Float3 operator*(const Float3& point, const Mat3& matrix)
	{
		return Float3{
			point.x * matrix.GetVal(0) + point.y * matrix.GetVal(1) + point.z * matrix.GetVal(2),
			point.x * matrix.GetVal(3) + point.y * matrix.GetVal(4) + point.z * matrix.GetVal(5),
			point.x * matrix.GetVal(6) + point.y * matrix.GetVal(7) + point.z * matrix.GetVal(8)
		};
	}
	constexpr void operator*=(Float3& point, const Mat3& matrix)
	{
		point.x = point.x * matrix.GetVal(0) + point.y * matrix.GetVal(1) + point.z * matrix.GetVal(2);
		point.y = point.x * matrix.GetVal(3) + point.y * matrix.GetVal(4) + point.z * matrix.GetVal(5);
		point.z = point.x * matrix.GetVal(6) + point.y * matrix.GetVal(7) + point.z * matrix.GetVal(8);
	}

	template<size_t columns, size_t rows>
	constexpr Matrix<rows, columns> Transpose(const Matrix<columns, rows>& matrix)
	{
		Matrix<rows, columns> retMat;
		for (size_t c = 0; c < columns; c++)
		{
			for (size_t r = 0; r < rows; r++)
				retMat.SetVal(r, c, matrix.GetVal(c, r));
		}
		return retMat;
	}
	template<size_t dim>
	constexpr float Determinant(const Matrix<dim>& matrix)
	{
		float num = 0;
		for (size_t mCol = 0; mCol < dim; mCol++)
		{
			Matrix<dim - 1> tempMat;
			size_t wCol = 0;

			for (size_t col = 0; col < dim; col++)
			{
				if (col == mCol)
				{
					continue;
				}
				for (size_t row = 0; row < dim - 1; row++)
				{
					tempMat.SetVal(row, wCol, matrix.GetVal(col, row + 1));
				}
				wCol++;
			}

			if (mCol % 2 == 0)
				num += matrix.GetVal(mCol, 0) * Determinant(tempMat);
			else
				num -= matrix.GetVal(mCol, 0) * Determinant(tempMat);
		}
		return num;
	}
	template<>
	constexpr float Determinant<2>(const Mat2& matrix)
	{
		return matrix.GetVal(0) * matrix.GetVal(3) - matrix.GetVal(1) * matrix.GetVal(2);
	}
	template<>
	constexpr float Determinant<3>(const Mat3& matrix)
	{
		float pos1 = matrix.GetVal(0, 0) * matrix.GetVal(1, 1) * matrix.GetVal(2, 2);
		float pos2 = matrix.GetVal(1, 0) * matrix.GetVal(2, 1) * matrix.GetVal(0, 2);
		float pos3 = matrix.GetVal(2, 0) * matrix.GetVal(0, 1) * matrix.GetVal(1, 2);

		float neg1 = matrix.GetVal(0, 2) * matrix.GetVal(1, 1) * matrix.GetVal(2, 0);
		float neg2 = matrix.GetVal(1, 2) * matrix.GetVal(2, 1) * matrix.GetVal(0, 0);
		float neg3 = matrix.GetVal(2, 2) * matrix.GetVal(0, 1) * matrix.GetVal(1, 0);

		return pos1 + pos2 + pos3 - neg1 - neg2 - neg3;
	}
	/*template<>
	constexpr float Determinant<4>(const Mat4& matrix)
	{
		float pos1 = matrix.GetVal(0, 0) * matrix.GetVal(1, 1) * matrix.GetVal(2, 2) * matrix.GetVal(3, 3);
		float pos2 = matrix.GetVal(1, 0) * matrix.GetVal(2, 1) * matrix.GetVal(3, 2) * matrix.GetVal(0, 3);
		float pos3 = matrix.GetVal(2, 0) * matrix.GetVal(3, 1) * matrix.GetVal(0, 2) * matrix.GetVal(1, 3);
		float pos4 = matrix.GetVal(3, 0) * matrix.GetVal(0, 1) * matrix.GetVal(1, 2) * matrix.GetVal(2, 3);

		float neg1 = matrix.GetVal(0, 3) * matrix.GetVal(1, 2) * matrix.GetVal(2, 1) * matrix.GetVal(3, 0);
		float neg2 = matrix.GetVal(1, 3) * matrix.GetVal(2, 2) * matrix.GetVal(3, 1) * matrix.GetVal(0, 0);
		float neg3 = matrix.GetVal(2, 3) * matrix.GetVal(3, 2) * matrix.GetVal(0, 1) * matrix.GetVal(1, 0);
		float neg4 = matrix.GetVal(3, 3) * matrix.GetVal(0, 2) * matrix.GetVal(1, 1) * matrix.GetVal(2, 0);

		return pos1 + pos2 + pos3 + pos4 - neg1 - neg2 - neg3 - neg4;
	}*/

	template<size_t dim>
	constexpr Matrix<dim> MatrixOfMinors(const Matrix<dim>& matrix)
	{
		Matrix<dim> retMat;
		for (size_t index = 0; index < dim * dim; index++)
		{
			size_t mColumn = floor(index / dim);
			size_t mRow = index % dim;

			Matrix<dim - 1> tempMat;

			size_t wRow = 0;

			for (size_t row = 0; row < dim; row++)
			{
				size_t wCol = 0;

				if (row == mRow)
					continue;

				for (size_t col = 0; col < dim; col++)
				{
					if (col == mColumn)
						continue;

					tempMat.SetVal(wCol, wRow, matrix.GetVal(col, row));
					wCol++;
				}
				wRow++;
			}

			retMat.SetVal(mColumn, mRow, Determinant(tempMat));
		}

		return retMat;
	}
	template<size_t dim>
	constexpr Matrix<dim> MatrixOfCofactors(const Matrix<dim>& matrix)
	{
		Matrix<dim> retMat = MatrixOfMinors(matrix);
		for (size_t i = 0; i < dim * dim; i++)
		{
			size_t column = floor(i / dim);
			size_t row = i % dim;

			retMat.SetVal(i, matrix.GetVal(i) * PowNeg1(row + column + 2));
		}
		return retMat;
	}

	template <size_t dim>
	constexpr Matrix<dim> Inverse(const Matrix<dim>& matrix)
	{
		float det = Determinant(matrix);
		Matrix<dim> retMat;

		if (det == 0)
			return retMat;

		retMat = MatrixOfMinors(matrix);
		retMat = MatrixOfCofactors(retMat);
		retMat = Transpose(retMat);
		retMat *= (1 / det);
		return retMat;
	}
	template<>
	constexpr Matrix<2> Inverse(const Matrix<2>& matrix)
	{
		float det = Determinant(matrix);
		float matArray[4] = {
			matrix.GetVal(3),
			-matrix.GetVal(1),
			-matrix.GetVal(2),
			matrix.GetVal(0)
		};
		Mat2 retMat;
		retMat.SetMatrix(matArray);

		retMat *= (1 / det);
		return retMat;
	}

	template <size_t cols, size_t rows>
	constexpr Matrix<cols, rows> GaussElim(const Matrix<cols, rows>& valMat, const Matrix<rows>& cMat)
	{
		Matrix<rows> calcMat{cMat};
		Matrix<cols, rows> eMat{valMat};

		// Get echelon form
		for (size_t i = 1; i < rows; i++)
		{
			for (size_t j = i; j < rows; j++)
			{
				Matrix<rows, 1> refRow = calcMat.GetRow(i - 1);
				Matrix<rows, 1> calcRow = calcMat.GetRow(j);
				float multiplier = calcRow.GetVal(i - 1) / refRow.GetVal(i - 1);
				refRow *= multiplier;
				calcRow = calcRow - refRow;

				Matrix<cols, 1> eRefRow = eMat.GetRow(i - 1) * multiplier;
				Matrix<cols, 1> eRow = eMat.GetRow(j);
				eRow = eRow - eRefRow;

				calcMat.SetRow(j, calcRow);
				eMat.SetRow(j, eRow);
			}
		}

		// Do the Rest
		for (size_t i = rows - 1; i >= 0 && i != SIZE_MAX; i--)
		{
			Matrix <rows, 1> rRow = calcMat.GetRow(i);
			Matrix <cols, 1> eRRow = eMat.GetRow(i);
			float mult = 1 / rRow.GetVal(i);
			rRow *= mult;
			eRRow *= mult;
			calcMat.SetRow(i, rRow);
			eMat.SetRow(i, eRRow);

			for (size_t j = i; j > 0; j--)
			{
				Matrix<rows, 1> refRow = calcMat.GetRow(i);
				Matrix<rows, 1> calcRow = calcMat.GetRow(j - 1);
				float multiplier = calcRow.GetVal(i);
				refRow *= multiplier;
				calcRow = calcRow - refRow;

				Matrix<cols, 1> eRefRow = eMat.GetRow(i) * multiplier;
				Matrix<cols, 1> eRow = eMat.GetRow(j - 1);
				eRow = eRow - eRefRow;

				calcMat.SetRow(j - 1, calcRow);
				eMat.SetRow(j - 1, eRow);
			}
		}

		return eMat;
	}

	template<size_t dim>
	constexpr Matrix<dim> GaussJordan(const Matrix<dim>& mat)
	{
		Matrix<dim> iMat{1};

		return GaussElim(iMat, mat);
	}

	constexpr Mat2 RotationMatrix2D(const float rotation)
	{
		Mat2 retMat;
		retMat.SetVal(0, cosf(Deg2Rad(rotation)));
		retMat.SetVal(1, -sinf(Deg2Rad(rotation)));
		retMat.SetVal(2, -retMat.GetVal(1));
		retMat.SetVal(3, retMat.GetVal(0));

		return retMat;
	}

	constexpr Mat2 Vec2Mat(const Float2& v1, const Float2& v2)
	{
		Mat2 mat;
		mat.SetVal(0, v1.x);
		mat.SetVal(1, v1.y);
		mat.SetVal(2, v2.x);
		mat.SetVal(3, v2.y);
		return mat;
	}
	constexpr Mat3 Vec2Mat(const Float3& v1, const Float3& v2, const Float3& v3)
	{
		float values[9] = {
			v2.x, v2.y, v2.z,
			v3.x, v3.y, v3.z,
			v1.x, v1.y, v1.z
		};

		Mat3 mat;
		mat.SetMatrix(values);
		return mat;
	}

	constexpr Mat4 MatPerspective(const float fovY, const float aspect, const float zNear, const float zFar)
	{
		Mat4 perspMat{0};

		float d = 1 / tan(MathLib::Deg2Rad(fovY / 2.0));
		perspMat[0][0] = d / aspect;
		perspMat[1][1] = d;
		perspMat[2][2] = (zNear + zFar) / (zNear - zFar);
		perspMat[2][3] = -1;
		perspMat[3][2] = 2 * zNear * zFar / (zNear - zFar);

		return perspMat;
	}
}
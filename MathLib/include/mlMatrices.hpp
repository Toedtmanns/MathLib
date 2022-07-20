#pragma once
#include "mlCommon.hpp"
#include <string.h>

namespace MathLib
{
	// Matrix maths

	template <size_t columns, size_t rows = columns>
	class EXPORT Matrix
	{
	protected:
		double m_Matrix[columns * rows]{0};

	public:
		constexpr Matrix(Matrix<columns, rows>&& other) noexcept
		{
			memcpy(m_Matrix, other.m_Matrix, sizeof(m_Matrix));
		}
		constexpr Matrix(const Matrix<columns, rows>& other)
		{
			memcpy(m_Matrix, other.m_Matrix, sizeof(m_Matrix));
		}
		constexpr Matrix()
		{

		}
		constexpr Matrix(const double defaultVal)
		{
			if (columns != rows)
				return;

			for (size_t i = 0; i < columns * rows; i += columns + 1)
			{
				m_Matrix[i] = defaultVal;
			}
		}

		constexpr void SetVal(const size_t column, const size_t row, const double value)
		{
			m_Matrix[column * rows + row] = value;
		}
		constexpr void SetVal(const size_t index, const double value)
		{
			m_Matrix[index] = value;
		}
		constexpr void SetMatrix(const double* matValues)
		{
			memcpy(m_Matrix, matValues, sizeof(m_Matrix));
		}

		// Getter methods

		constexpr const double GetVal(const size_t column, const size_t row) const
		{
			return m_Matrix[column * rows + row];
		}
		constexpr const double GetVal(const size_t index) const
		{
			return m_Matrix[index];
		}
		constexpr double* GetMatrix()
		{
			return m_Matrix;
		}
		constexpr const double* GetMatrix() const
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
		double* operator[](const size_t index) const
		{
			double retArr[rows];
			memcpy(retArr, m_Matrix + index * rows, rows);
			return retArr;
		}

		constexpr Matrix<columns, rows> operator*(const double value)
		{
			Matrix<columns, rows> retMat{0};
			for (size_t i = 0; i < columns * rows; i++)
			{
				retMat.SetVal(i, m_Matrix[i] * value);
			}
			return retMat;
		}
		template<size_t oColumns>
		constexpr Matrix<oColumns, rows> operator*(const Matrix<oColumns, columns>& other) const
		{
			Matrix<oColumns, rows> retMat;
			for (size_t c = 0; c < oColumns; c++)
			{
				for (size_t r = 0; r < rows; r++)
				{
					double num = 0;
					for (size_t calcCol = 0; calcCol < columns; calcCol++)
					{
						num = num + m_Matrix[calcCol * rows + r] * other.GetMatrix()[calcCol * columns + c];
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
					double num = 0;
					for (size_t calcCol = 0; calcCol < columns; calcCol++)
					{
						num = num + m_Matrix[calcCol * rows + r] * other.GetMatrix()[calcCol * columns + c];
					}
					m_Matrix[c * rows + r] = num;
				}
			}
		}
		constexpr Matrix<columns, rows> operator+(const double value) const
		{
			Matrix<columns, rows> retMat;
			for (size_t i = 0; i < columns * rows; i++)
				retMat.SetVal(i, m_Matrix[i] + value);
			return retMat;
		}
		constexpr void operator+=(const double value)
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
		constexpr Matrix<columns, rows> operator-(const double value) const
		{
			Matrix<columns, rows> retMat;
			for (size_t i = 0; i < columns * rows; i++)
				retMat.SetVal(i, m_Matrix[i] - value);
			return retMat;
		}
		constexpr void operator-=(const double value)
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

	EXPORT typedef Matrix<4> Mat4;
	EXPORT typedef Matrix<3> Mat3;
	EXPORT typedef Matrix<2> Mat2;
	EXPORT typedef Matrix<1, 2> Mat1x2;
	EXPORT typedef Matrix<1, 3> Mat1x3;
	EXPORT typedef Matrix<1, 4> Mat1x4;
	EXPORT typedef Matrix<2, 1> Mat2x1;
	EXPORT typedef Matrix<2, 3> Mat2x3;
	EXPORT typedef Matrix<2, 4> Mat2x4;
	EXPORT typedef Matrix<3, 1> Mat3x1;
	EXPORT typedef Matrix<3, 2> Mat3x2;
	EXPORT typedef Matrix<3, 4> Mat3x4;
	EXPORT typedef Matrix<4, 1> Mat4x1;
	EXPORT typedef Matrix<4, 2> Mat4x2;
	EXPORT typedef Matrix<4, 3> Mat4x3;

	template<size_t columns, size_t rows>
	EXPORT constexpr Matrix<rows, columns> Transpose(const Matrix<columns, rows>& matrix)
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
	EXPORT constexpr double Determinant(const Matrix<dim>& matrix)
	{
		double num = 0;
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
	EXPORT constexpr double Determinant<2>(const Mat2& matrix)
	{
		return matrix.GetVal(0) * matrix.GetVal(3) - matrix.GetVal(1) * matrix.GetVal(2);
	}
	template<>
	EXPORT constexpr double Determinant<3>(const Mat3& matrix)
	{
		double pos1 = matrix.GetVal(0, 0) * matrix.GetVal(1, 1) * matrix.GetVal(2, 2);
		double pos2 = matrix.GetVal(1, 0) * matrix.GetVal(2, 1) * matrix.GetVal(0, 2);
		double pos3 = matrix.GetVal(2, 0) * matrix.GetVal(0, 1) * matrix.GetVal(1, 2);

		double neg1 = matrix.GetVal(0, 2) * matrix.GetVal(1, 1) * matrix.GetVal(2, 0);
		double neg2 = matrix.GetVal(1, 2) * matrix.GetVal(2, 1) * matrix.GetVal(0, 0);
		double neg3 = matrix.GetVal(2, 2) * matrix.GetVal(0, 1) * matrix.GetVal(1, 0);

		return pos1 + pos2 + pos3 - neg1 - neg2 - neg3;
	}

	template<size_t dim>
	EXPORT constexpr Matrix<dim> MatrixOfMinors(const Matrix<dim>& matrix)
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
	}

	template<size_t dim>
	EXPORT constexpr Matrix<dim> MatrixOfCofactors(const Matrix<dim>& matrix)
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
	EXPORT constexpr Matrix<dim> Inverse(const Matrix<dim>& matrix)
	{
		double det = Determinant(matrix);
		Matrix<dim> retMat;

		if (det == 0)
			return retMat;

		retMat = MatrixOfMinors(matrix);
		retMat = MatrixOfCofactors(retMat);
		retMat = Transpose(retMat);
		retMat = retMat * (1 / det);
		return retMat;
	}
	template<>
	EXPORT constexpr Matrix<2> Inverse(const Matrix<2>& matrix)
	{
		double matArray[4] = {
			matrix.GetVal(3),
			-matrix.GetVal(1),
			-matrix.GetVal(2),
			matrix.GetVal(0)
		};
		Mat2 retMat;
		retMat.SetMatrix(matArray);
		return retMat;
	}

	EXPORT constexpr Mat2 RotationMatrix2D(const double rotation)
	{
		Mat2 retMat;
		retMat.SetVal(0, cos(Deg2Rad(rotation)));
		retMat.SetVal(1, -sin(Deg2Rad(rotation)));
		retMat.SetVal(2, -retMat.GetVal(1));
		retMat.SetVal(3, retMat.GetVal(1));

		return retMat;
	}

	// Old matrix classes

	
	class DEPRECATEDCLASS("Integer matrices are unnecessary and reworked matrices are available") EXPORT MatrixI
	{
		unsigned int m_Rows, m_Columns;
		int** m_Matrix;

	public:
		MatrixI(MatrixI&& other) noexcept;
		MatrixI(const MatrixI& other);
		MatrixI(const unsigned int& dim);
		MatrixI(const unsigned int& columns, const unsigned int& rows);
		MatrixI(const unsigned int& columns, const unsigned int& rows, const int** const columnArray);

		void SetNum(const unsigned int& column, const unsigned int& row, const int& value);
		void SetColumn(const unsigned int& column, const int* const content);
		void SetMatrix(const int** const matrix);

		const int& GetNum(const unsigned int& column, const unsigned int& row) const;
		int* GetColumn(const unsigned int& column) const;
		int* GetArray() const;
		int** GetMatrix() const;
		const unsigned int& GetRowCount() const;
		const unsigned int& GetColumnCount() const;

		void operator=(MatrixI&& other) noexcept;
		void operator=(const MatrixI& other);
		MatrixI operator+(const int& value) const;
		MatrixI operator-(const int& value) const;
		MatrixI operator*(const int& value) const;
		MatrixI operator/(const int& value) const;
		MatrixI operator+(const MatrixI& other) const;
		MatrixI operator-(const MatrixI& other) const;
		MatrixI operator*(const MatrixI& other) const;
		int* operator[](const unsigned int& column);
		const int* operator[](const unsigned int& column) const;

		~MatrixI();
	};

	class EXPORT MatrixF
	{
		double** m_Matrix;
		unsigned int m_Rows, m_Columns;

	public:
		MatrixF(MatrixF&& other) noexcept;
		MatrixF(const MatrixF& other);
		MatrixF(const unsigned int& dim);
		MatrixF(const unsigned int& columns, const unsigned int& rows);
		MatrixF(const unsigned int& columns, const unsigned int& rows, const double** const columnArray);

		void SetNum(const unsigned int& column, const unsigned int& row, const double& value);
		void SetColumn(const unsigned int& column, const double* const content);
		void SetMatrix(const double** const matrix);

		const double& GetNum(const unsigned int& column, const unsigned int& row) const;
		double* GetColumn(const unsigned int& column) const;
		double* GetArray() const;
		double** GetMatrix() const;
		const unsigned int& GetRowCount() const;
		const unsigned int& GetColumnCount() const;

		void operator=(MatrixF&& other) noexcept;
		void operator=(const MatrixF& other);
		MatrixF operator+(const double& value) const;
		MatrixF operator-(const double& value) const;
		MatrixF operator*(const double& value) const;
		MatrixF operator/(const double& value) const;
		MatrixF operator+(const MatrixF& other) const;
		MatrixF operator-(const MatrixF& other) const;
		MatrixF operator*(const MatrixF& other) const;
		double* operator[](const unsigned int& column);
		const double* operator[](const unsigned int& column) const;

		~MatrixF();
	};

	DEPRECATED()
	EXPORT void PrintContent(const MatrixI& mat);
	DEPRECATED()
	EXPORT void PrintProperties(const MatrixI& mat);
	DEPRECATED()
	EXPORT int MatrixGetDet(const MatrixI& mat);
	DEPRECATED()
	EXPORT MatrixI MatrixOfMinors(const MatrixI& mat);
	DEPRECATED()
	EXPORT MatrixI MatrixOfCofactors(const MatrixI& mat);
	DEPRECATED()
	EXPORT MatrixI MatrixTranspose(const MatrixI& mat);
	DEPRECATED()
	EXPORT MatrixI MatrixAdjugate(const MatrixI& mat);
	DEPRECATED()
	EXPORT MatrixF MatrixInverse(const MatrixI& mat);

	EXPORT void PrintContent(const MatrixF& mat);
	EXPORT void PrintProperties(const MatrixF& mat);
	EXPORT double MatrixGetDet(const MatrixF& mat);
	EXPORT MatrixF MatrixOfMinors(const MatrixF& mat);
	EXPORT MatrixF MatrixOfCofactors(const MatrixF& mat);
	EXPORT MatrixF MatrixTranspose(const MatrixF& mat);
	EXPORT MatrixF MatrixAdjugate(const MatrixF& mat);
	EXPORT MatrixF MatrixInverse(const MatrixF& mat);
}